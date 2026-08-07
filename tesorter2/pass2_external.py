"""
pass2_external.py — helpers for the --pass2-classified-fasta feature.

Ported from github.com/cwb14/TEsorter branch `my-new-idea2` (TEsorter/app.py,
head commit b398509). The upstream helpers used Biopython SeqIO and TEsorter's
CommonClassification namedtuple. Here we use pyfastx (already a TEsorter2
dependency) and emit dicts matching TEsorter2's classifications-dict shape
(`id/order/superfamily/clade/complete/strand/domains/score/secondary`).
"""

import logging
import os
import re

import pyfastx

log = logging.getLogger(__name__)


_COORD_HEADER_RE = re.compile(
    r'^(?P<id>\S+?:\d+[-\.]+\d+)#(?P<order>[^/]+)/(?P<sfam>[^/]+)/(?P<clade>\S+)$'
)


def _format_gff_id(s):
    """TEsorter's format_gff_id: strip anything after a '#'. Trivial but kept
    as a named helper so the intent reads."""
    return s.split("#", 1)[0]


def parse_cls_from_fasta_header(header):
    """Parse `>id#Order/Superfamily/Clade` -> (id, order, superfamily, clade).

    Returns None for headers without a '#' or without at least Order/Superfamily.
    Missing slots are filled with 'Unknown'/'unknown'.
    """
    h = header.strip()
    if h.startswith('>'):
        h = h[1:]
    h = h.split(None, 1)[0]

    if '#' not in h:
        return None
    raw_id, cls = h.split('#', 1)
    raw_id = _format_gff_id(raw_id)

    parts = cls.split('/')
    if len(parts) < 2:
        return None

    order = parts[0] or 'Unknown'
    superfamily = parts[1] or 'unknown'
    clade = parts[2] if len(parts) >= 3 and parts[2] else 'unknown'
    return raw_id, order, superfamily, clade


def extend_hmm_classifications_from_fasta(hmm_cls, fasta_path, db_seq_to_dbs):
    """Merge external FASTA classifications into hmm_cls in place.

    Each added entry is shaped like TEsorter2's reconciled classifications
    dict so classifier.py:568 still sees the fields it expects. Also extends
    db_seq_to_dbs so classify_from_blast accepts hits pointing at these IDs.
    """
    if fasta_path is None:
        return

    added = 0
    skipped = 0
    # pyfastx is faster than Biopython and already in TEsorter2's deps.
    fa = pyfastx.Fasta(fasta_path, build_index=True, full_name=True)
    for rec in fa:
        parsed = parse_cls_from_fasta_header(rec.name)
        if not parsed:
            skipped += 1
            continue
        sid, order, superfamily, clade = parsed
        if sid in hmm_cls:
            continue
        hmm_cls[sid] = {
            "id": sid,
            "order": order,
            "superfamily": superfamily,
            "clade": clade,
            "complete": "none",
            "strand": "?",
            "domains": "none",
            "score": 0.0,
            "secondary": [],
        }
        db_seq_to_dbs.setdefault(sid, set()).add("external")
        added += 1

    log.info(
        f"extended pass-1 classifications with {added} entries from {fasta_path} "
        f"({skipped} headers skipped: not parseable)"
    )


def merge_classified_fastas(out_fa, fa_primary, fa_extra=None, clean_nucl=True):
    """Write out_fa as the combined pass-2 target FASTA.

    IDs are stripped of any trailing '#...' annotation so they round-trip
    cleanly through mmseqs. Primary takes precedence on duplicate IDs.
    When clean_nucl=True, non-ATCG characters are stripped before writing.
    """
    seen = set()
    n = 0

    def emit_records(path, fout):
        nonlocal n
        fa = pyfastx.Fasta(path, build_index=True)
        for rec in fa:
            rid = _format_gff_id(rec.name)
            if rid in seen:
                continue
            seen.add(rid)
            seq = str(rec.seq)
            if clean_nucl:
                seq = "".join(c for c in seq.upper() if c in "ATCG")
            fout.write(f">{rid}\n{seq}\n")
            n += 1

    with open(out_fa, "w") as fout:
        emit_records(fa_primary, fout)
        if fa_extra:
            emit_records(fa_extra, fout)

    suffix = " [non-ATCG stripped]" if clean_nucl else ""
    log.info(f"pass-2 database FASTA written: {out_fa} ({n} unique IDs){suffix}")


def update_classified_fasta_headers(fasta_path, hmm_cls, tmpdir):
    """Rewrite a pass2-classified FASTA, upgrading `unknown` slots from hmm_cls.

    Only headers shaped `>chr:start-end#Order/Superfamily/Clade` where the ID
    matches an entry in hmm_cls AND at least one of Order/Superfamily/Clade
    contains 'unknown' are candidates. Returns the path to the written FASTA,
    or None if the input was None.
    """
    if fasta_path is None:
        return None

    os.makedirs(tmpdir, exist_ok=True)
    updated_path = os.path.join(tmpdir, "pass2_classified_updated.fa")
    n_updated = 0
    n_total = 0

    fa = pyfastx.Fasta(fasta_path, build_index=True, full_name=True)
    with open(updated_path, "w") as fout:
        for rec in fa:
            n_total += 1
            header = rec.name.split(None, 1)[0]
            seq = str(rec.seq)

            if '#' in header:
                raw_id_part, _ = header.split('#', 1)
                rid = _format_gff_id(raw_id_part)

                m = _COORD_HEADER_RE.match(header)
                if m and rid in hmm_cls:
                    old_order = m.group('order')
                    old_sfam = m.group('sfam')
                    old_clade = m.group('clade')

                    if 'unknown' in (old_order.lower(), old_sfam.lower(), old_clade.lower()):
                        cls = hmm_cls[rid]
                        new_cls = "{}/{}/{}".format(
                            cls["order"], cls["superfamily"], cls["clade"]
                        )
                        header = f"{raw_id_part}#{new_cls}"
                        n_updated += 1

            fout.write(f">{header}\n{seq}\n")

    log.info(
        f"updated {n_updated}/{n_total} headers in pass2-classified-fasta "
        f"using pass-1 classifications"
    )
    return updated_path


def clean_fasta_atcg(path):
    """In-place ATCG-only cleaner. Kept for completeness / upstream parity.

    Not actively called by TEsorter2's pass-2 because the mmseqs wrapper
    already cleans the query and merge_classified_fastas cleans the DB.
    """
    tmp = path + ".atcg_clean.tmp"
    with open(path) as fin, open(tmp, "w") as fout:
        buf = []
        header = [None]

        def flush():
            if header[0] is None:
                return
            seq = "".join(buf).upper()
            seq = "".join(c for c in seq if c in "ATCG")
            fout.write(header[0] + "\n")
            fout.write(seq + "\n")

        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header[0] = line
                buf.clear()
            else:
                buf.append(line)
        flush()
    os.replace(tmp, path)
    log.info(f"non-ATCG characters removed from {path}")
