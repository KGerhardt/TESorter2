"""
tesorter_output.py — Generate TEsorter-compatible output files.

Produces the same output files as TEsorter from our SQLite database:
  - {prefix}.cls.tsv  — TE classifications (already produced by classifier)
  - {prefix}.dom.gff3 — Domain annotations in GFF3 format
  - {prefix}.dom.tsv  — Domain hit summary table
  - {prefix}.dom.faa  — Domain protein sequences in FASTA
  - {prefix}.cls.lib  — RepeatMasker library (classified sequences)
  - {prefix}.cls.pep  — Classified protein domain sequences
"""

import logging
import re

import pyfastx

from .sequence import (parse_frame_suffix, aa_to_nucl_coords,
                      load_sequences_dict, open_input, revcomp)
from .classifier import parse_clade_rexdb, parse_clade_gydb

log = logging.getLogger(__name__)


def _parse_gene_clade(model_name, db_name):
    """Extract gene and clade from model name."""
    if db_name.startswith("rexdb") or db_name in ("line", "tir"):
        if ":" in model_name:
            gene = model_name.split(":")[1]
            clade = model_name.split("/")[-1].split(":")[0]
        else:
            gene = model_name
            clade = model_name
    elif db_name == "gydb":
        parts = model_name.split("_", 1)
        gene = parts[0]
        clade = parts[1] if len(parts) > 1 else model_name
    elif db_name.startswith("sine"):
        gene = "SINE"
        clade = "SINE"
    else:
        gene = model_name
        clade = model_name
    return gene, clade


def _format_gff_id(seq_id):
    """Make ID safe for GFF3 (replace special chars)."""
    return re.sub(r"[;=|]", "_", seq_id)


# Hit tables these exports can read. Table names cannot be bound parameters,
# so the name is validated against this allowlist before being interpolated.
_HITS_TABLES = ("legacy_hits",)


def _hits_query(columns, hits_table):
    """Build the per-database hit query used by every domain-level export."""
    if hits_table not in _HITS_TABLES:
        raise ValueError(f"unknown hits table: {hits_table}")
    return f"""
        SELECT {columns}
        FROM {hits_table}
        WHERE database = ?
        ORDER BY base_seq, env_from
    """


def domain_keys(hits, indices):
    """Row keys for the domain assignments the classifier selected.

    Pass the result as `keep` to the domain exports so they emit one row per
    region -- the model that won hmm2best -- instead of every model that hit
    it. Without this they fall back to threshold filtering alone, which leaves
    all the losing models in.
    """
    return {(str(hits["target"][i]), str(hits["model"][i]),
             int(hits["env_from"][i]), int(hits["env_to"][i]))
            for i in indices}


def _fmt_cls(*levels):
    """Join classification levels the way TEsorter's fmt_cls does.

    Drops 'unknown' and any level that repeats one already emitted; 'mixture'
    is kept, so a conflicted call comes out as e.g. 'LTR/Gypsy/mixture' or
    plain 'mixture' when the conflict reaches the order level.
    """
    values = []
    for level in levels:
        if level == "unknown" or level in values:
            continue
        values.append(level)
    return "/".join(values)


def _write_fasta(f, header, seq, wrap=None):
    """Write one FASTA record, wrapping the sequence if TEsorter wraps it."""
    f.write(f">{header}\n")
    if not wrap:
        f.write(f"{seq}\n")
        return
    for i in range(0, len(seq), wrap):
        f.write(f"{seq[i:i + wrap]}\n")


def _skip_row(keep, target, model, env_from, env_to,
              hmm_cov, evalue, acc, norm_score):
    """Whether a raw hit row should be left out of a domain export."""
    if keep is not None:
        return (target, model, env_from, env_to) not in keep
    return not (hmm_cov >= 20 and evalue <= 1e-3 and acc >= 0.5
                and norm_score >= 0.1)


def _domain_rows(conn, db_name, nucl_lengths, seq_type="nucl",
                 hits_table="legacy_hits", keep=None):
    """Filtered domain hits, in the order TEsorter emits them.

    TEsorter writes every domain-level file from one list sorted by sequence
    then genomic start, so a minus-strand element's domains come out in
    ascending coordinate order -- the reverse of the translated-frame env_from
    order the hits are stored in. Resolving coordinates once here keeps the
    filtering and the ordering identical across .dom.gff3, .dom.tsv, .dom.faa
    and .cls.pep.
    """
    rows = conn.execute(_hits_query(
        """target_name, query_name, dom_score, i_evalue, acc,
           hmm_from, hmm_to, query_len, env_from, env_to,
           base_seq, strand, frame, domain_type""",
        hits_table), (db_name,)).fetchall()

    out = []
    for (target, model, score, evalue, acc, hmm_from, hmm_to, model_len,
         env_from, env_to, base_seq, strand, frame_num, domain_type) in rows:
        hmm_cov = 100.0 * (hmm_to - hmm_from + 1) / model_len
        norm_score = round(score / model_len, 2)

        if _skip_row(keep, target, model, env_from, env_to,
                     hmm_cov, evalue, acc, norm_score):
            continue

        # Convert to nucleotide coordinates
        if seq_type == "nucl" and strand in ("+", "-"):
            nuc_start, nuc_end = aa_to_nucl_coords(
                env_from, env_to, strand, int(frame_num),
                nucl_lengths.get(base_seq, 0))
            nuc_frame = int(frame_num)
        else:
            nuc_start, nuc_end = env_from, env_to
            nuc_frame = "."

        gene, clade = _parse_gene_clade(model, db_name)
        gid = f"{_format_gff_id(base_seq)}|{model}"
        out.append({
            "target": target, "base_seq": base_seq, "gid": gid,
            "env_from": env_from, "env_to": env_to,
            "nuc_start": nuc_start, "nuc_end": nuc_end,
            "strand": strand, "frame": nuc_frame,
            "norm_score": norm_score, "evalue": evalue, "acc": acc,
            "hmm_cov": hmm_cov, "clade": clade, "domain_type": domain_type,
            "name": f"{clade}-{gene.split('-')[-1] if '-' in gene else gene}",
        })

    out.sort(key=lambda r: (r["base_seq"], r["nuc_start"]))
    return out


def _domain_attr(row):
    """The GFF3 attribute string shared by .dom.gff3, .dom.faa and .cls.pep."""
    return (f"ID={row['gid']};Name={row['name']};gene={row['domain_type']};"
            f"clade={row['clade']};coverage={row['hmm_cov']:.1f};"
            f"evalue={row['evalue']};probability={row['acc']}")


def generate_dom_gff3(conn, prefix, db_name, nucl_lengths, seq_type="nucl",
                      hits_table="legacy_hits", keep=None, rows=None):
    """Generate domain annotation GFF3 file.

    Format: chr  TEsorter  CDS  start  end  score  strand  frame  attributes
    """
    out_path = f"{prefix}.dom.gff3"

    if rows is None:
        rows = _domain_rows(conn, db_name, nucl_lengths, seq_type,
                            hits_table, keep)

    with open(out_path, "w") as f:
        for row in rows:
            line = [row["base_seq"], "TEsorter", "CDS",
                    str(row["nuc_start"]), str(row["nuc_end"]),
                    f"{row['norm_score']}", row["strand"], str(row["frame"]),
                    _domain_attr(row)]
            f.write("\t".join(line) + "\n")

    log.info(f"  {out_path}: {len(rows)} domain annotations")
    return out_path


def generate_dom_tsv(conn, prefix, db_name, nucl_lengths=None,
                     seq_type="nucl", hits_table="legacy_hits", keep=None,
                     rows=None):
    """Generate domain hit summary table."""
    out_path = f"{prefix}.dom.tsv"

    if rows is None:
        rows = _domain_rows(conn, db_name, nucl_lengths or {}, seq_type,
                            hits_table, keep)

    with open(out_path, "w") as f:
        f.write("#id\tlength\tevalue\tcoverge\tprobability\tscore\n")
        for row in rows:
            seq_len = row["env_to"] - row["env_from"] + 1
            f.write(f"{row['gid']}\t{seq_len}\t{row['evalue']}\t"
                    f"{row['hmm_cov']:.1f}\t{row['acc']}\t"
                    f"{row['norm_score']}\n")

    log.info(f"  {out_path}")
    return out_path


def generate_dom_faa(conn, prefix, db_name, aa_fasta, nucl_lengths=None,
                     seq_type="nucl", hits_table="legacy_hits", keep=None,
                     rows=None):
    """Generate domain protein sequences FASTA."""
    out_path = f"{prefix}.dom.faa"

    if rows is None:
        rows = _domain_rows(conn, db_name, nucl_lengths or {}, seq_type,
                            hits_table, keep)

    # Load translated sequences
    seqs = load_sequences_dict(aa_fasta)

    with open(out_path, "w") as f:
        for row in rows:
            if row["target"] not in seqs:
                continue
            subseq = seqs[row["target"]][row["env_from"] - 1:row["env_to"]]
            f.write(f">{row['gid']} {_domain_attr(row)}\n{subseq}\n")

    log.info(f"  {out_path}")
    return out_path


def generate_cls_lib(input_fasta, prefix, classifications,
                     no_reverse=False):
    """Generate RepeatMasker library FASTA.

    Every sequence gets its classification appended to the ID, keeping the
    original header as the description:
        >seqid#Order/Superfamily/Clade original_header
    Unclassified sequences are labelled '#Unknown'. Minus-strand sequences are
    reverse complemented unless no_reverse.
    """
    out_path = f"{prefix}.cls.lib"

    cls_dict = {r["id"]: r for r in classifications}

    # Load all sequences into memory
    fa = pyfastx.Fasta(input_fasta, build_index=True)
    all_seqs = {rec.name: str(rec.seq) for rec in fa}

    with open(out_path, "w") as f:
        for name, seq in all_seqs.items():
            if name in cls_dict:
                cl = cls_dict[name]
                strand = cl.get("strand", "+")
                cls_str = _fmt_cls(cl["order"], cl["superfamily"], cl["clade"])

                if not no_reverse and strand == "-":
                    seq = revcomp(seq)
            else:
                cls_str = "Unknown"

            base_id = name.split("#")[0]
            _write_fasta(f, f"{base_id}#{cls_str} {name}", seq, wrap=60)

    log.info(f"  {out_path}")
    return out_path


def generate_cls_pep(conn, prefix, db_name, aa_fasta, classifications,
                     nucl_lengths=None, seq_type="nucl",
                     hits_table="legacy_hits", keep=None, rows=None):
    """Generate classified protein domain sequences.

    Each domain sequence gets the full classification in the ID:
        >seqid#Order/Superfamily/Clade#gene|clade
    """
    out_path = f"{prefix}.cls.pep"

    if rows is None:
        rows = _domain_rows(conn, db_name, nucl_lengths or {}, seq_type,
                            hits_table, keep)

    cls_dict = {r["id"]: r for r in classifications}
    seqs = load_sequences_dict(aa_fasta)

    with open(out_path, "w") as f:
        for row in rows:
            base_seq = row["base_seq"]
            if row["target"] not in seqs or base_seq not in cls_dict:
                continue

            subseq = seqs[row["target"]][row["env_from"] - 1:row["env_to"]]
            cl = cls_dict[base_seq]
            clade, domain_type = row["clade"], row["domain_type"]

            cls_str = _fmt_cls(cl["order"], cl["superfamily"], cl["clade"])

            raw_id = base_seq.split("#")[0]
            new_id = f"{raw_id}#{cls_str}#{domain_type}|{clade}"

            gid = row["gid"]
            attr = (f"ID={gid};Name={clade}-{domain_type};"
                    f"gene={domain_type};clade={clade};"
                    f"coverage={row['hmm_cov']:.1f};evalue={row['evalue']};"
                    f"probability={row['acc']}")

            # TEsorter repeats the domain id as the description before the
            # attribute string, and wraps the peptide at 60 columns.
            _write_fasta(f, f"{new_id} {gid} {attr}", subseq, wrap=60)

    log.info(f"  {out_path}")
    return out_path


def generate_all_outputs(conn, prefix, db_name, input_fasta, aa_fasta,
                         nucl_lengths, classifications,
                         seq_type="nucl", no_reverse=False, no_library=False,
                         hits_table="legacy_hits", domain_files=True, keep=None):
    """Generate all TEsorter-compatible output files.

    domain_files=False emits only the library (.cls.lib), which needs nothing
    but the input sequences and their classifications. The domain-level files
    carry six-frame translated coordinates, so they are only meaningful when
    the run actually took the six-frame protein path -- not under --bath, where
    bathsearch reports nucleotide coordinates directly, and not for
    DNA-alphabet databases searched by nhmmer.
    """
    log.info(f"Generating TEsorter-format outputs for {db_name}")

    if domain_files:
        # One filtered, coordinate-sorted row list feeds every domain file, so
        # they stay consistent with each other and with TEsorter's ordering.
        rows = _domain_rows(conn, db_name, nucl_lengths, seq_type,
                            hits_table, keep)
        generate_dom_gff3(conn, prefix, db_name, nucl_lengths, seq_type,
                          rows=rows)
        generate_dom_tsv(conn, prefix, db_name, rows=rows)

        if aa_fasta:
            generate_dom_faa(conn, prefix, db_name, aa_fasta, rows=rows)
            generate_cls_pep(conn, prefix, db_name, aa_fasta, classifications,
                             rows=rows)

    if not no_library:
        generate_cls_lib(input_fasta, prefix, classifications, no_reverse)
