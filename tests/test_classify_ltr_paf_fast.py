import pytest

from tesorter2.classify_ltr_paf_fast import process_paf


def _paf(pieces, qlen, tlen, dvs=None):
    """PAF lines for query q against target t, one per (strand, qs, qe, ts, te)."""
    return ["q\t%d\t%d\t%d\t%s\tt\t%d\t%d\t%d\t%d\t%d\t0\ttp:A:P\tdv:f:%s\n"
            % (qlen, qs, qe, strand, tlen, ts, te, qe - qs, qe - qs, dv)
            for (strand, qs, qe, ts, te), dv in zip(pieces, dvs or [0.02] * len(pieces))]


@pytest.mark.parametrize("pieces, qlen, tlen, qcov, tcov, verdict", [
    pytest.param([("+", 15, 516, 398, 938), ("+", 4, 516, 1, 551), ("+", 4, 436, 774, 1243)],
                 517, 1251, 512 / 517, 550 / 1251, "fail", id="one-stretch-on-three-target-parts"),
    pytest.param([("+", 0, 1000, 0, 1000), ("+", 500, 1500, 1000, 2000)],
                 1500, 2000, 1000 / 1500, 1000 / 2000, "fail", id="in-order-but-reusing-500bp"),
    pytest.param([("+", 0, 1000, 0, 1000), ("+", 997, 2000, 1100, 2100)],
                 2000, 2100, 2000 / 2000, 2000 / 2100, "pass", id="broken-match-sharing-3bp"),
    pytest.param([("+", 0, 1000, 1000, 2000), ("+", 1000, 2000, 0, 1000)],
                 2000, 2000, 1000 / 2000, 1000 / 2000, "fail", id="out-of-order-on-target"),
    pytest.param([("-", 0, 1000, 1000, 2000), ("-", 1000, 2000, 0, 1000)],
                 2000, 2000, 2000 / 2000, 2000 / 2000, "pass", id="reverse-strand-runs-backwards"),
    pytest.param([("+", 0, 1000, 0, 1000), ("-", 1000, 2000, 1000, 2000)],
                 2000, 2000, 1000 / 2000, 1000 / 2000, "fail", id="opposite-strands"),
])
def test_coverage_adds_up_in_order_pieces_only(pieces, qlen, tlen, qcov, tcov, verdict):
    (row,) = process_paf(_paf(pieces, qlen, tlen), min_pid=0.7, min_qcov=0.7, min_tcov=0.7)
    assert row[3] == pytest.approx(qcov)
    assert row[4] == pytest.approx(tcov)
    assert row[1] == verdict


def test_identity_comes_from_the_chain_alone():
    """The chain is 70% identical over 300 bp, then 84% over 600 bp. A third
    alignment laid over most of the same query stretch, off the chain, must not
    lift the chain's 79.3% identity past 80%."""
    pieces = [("+", 0, 300, 0, 300), ("+", 300, 900, 300, 900), ("+", 100, 1000, 100, 700)]
    (row,) = process_paf(_paf(pieces, 1000, 1000, dvs=[0.3, 0.16, 0.0]),
                         min_pid=0.8, min_qcov=0.8, min_tcov=0.8)
    assert row[2:5] == pytest.approx(((0.7 * 300 + 0.84 * 600) / 900, 0.9, 0.9))
    assert row[1] == "fail"
