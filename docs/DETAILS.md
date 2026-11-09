# TEsorter2: technical details

The practical guide is the [README](../README.md). This page explains how TEsorter2 works.

## Databases

The HMM databases (REXdb, GyDB2, LINE, TIR, AnnoSINE, a Pfam-derived TE set and four Dfam curated
collections) ship inside the package, as they do in TEsorter, so there is no download step and no
configuration: `TEsorter2 input.fasta` works straight after install.

The Dfam collections ship gzipped. An HMM file is text and gives back ~84% of its size, and nothing
has to unpack them: pyHMMER's `HMMFile` detects gzip, and HMMER's `nhmmer` reads a `.gz` path
directly. So there is no first-use step, nothing is written into the package directory at run time,
and concurrent jobs sharing one install cannot race over a half-written file.

> **GyDB2 is shipped modified.** Its 2009-era profiles carried no `COMPO` record, which HMMER and
> pyHMMER tolerate but stricter HMMER3 readers (notably `nail`) reject outright. The bundled copy
> adds `COMPO` to all 314 profiles, derived from each model's own match emissions by HMMER's
> `p7_hmm_SetComposition()` — no probability was altered and no model was rebuilt. Verified to give
> byte-identical hits and classifications. Full notice, rationale and reproduction steps in
> [`tesorter2/database/LICENSES.md`](tesorter2/database/LICENSES.md).

To use a custom collection of HMM databases instead, point TEsorter2 at its directory:

```bash
TEsorter2 input.fasta --db-dir /path/to/db   # or: export TESORTER2_DB=/path/to/db
```

Individual databases can also be passed by path: `-db /path/to/custom.hmm` or, as in TEsorter,
`--db-hmm /path/to/custom.hmm`. A custom file is read with the classifier of the bundled database
whose name its file name contains (`my_gydb.hmm` reads as GyDB), otherwise with REXdb's naming
convention.

## Two modes

What the input *is* decides nearly everything downstream: which engines can run, and whether there
is a per-element classification at all. As in TEsorter, element mode is the default and `-genome`
selects genome mode:

```bash
TEsorter2 lib.fa                # pre-extracted TE sequences -> per-element classifications
TEsorter2 asm.fa -genome        # whole assembly -> domain-level GFF3 + summary
```


## Engines and the cascade

Amino-acid databases are searched by a staged cascade chosen with `--stages` (default
`nail,hmmer,bath`). Each stage searches only the sequences that earlier stages left unresolved, per
database, so the fast engine strips the easy majority and the sensitive ones see a shrinking
remainder. Naming a single engine (`--stages hmmer`, `--stages bath`) runs that engine alone. When
`--stages` is not given, a stage whose engine is not installed is skipped with a warning; an explicit
`--stages` list fails instead. DNA databases are always searched with nhmmer.

- **nail** ([TravisWheelerLab/nail](https://github.com/TravisWheelerLab/nail)) seeds with MMseqs2 and
  approximates HMMER's Forward/Backward. It reads HMMER3 `.hmm` files directly and is amino-acid only.
  Verified against nail 0.7.1 with MMseqs2 18.8cc5c. On many clusters `/usr/bin/nail` is the `mailx`
  mail client; TEsorter2 checks the binary and ignores it.
- **hmmer**: pyHMMER `hmmsearch` over six-frame translations, in-process.
- **bath**: `bathsearch --fs` on the nucleotide sequence (see below).

| | element mode (pre-extracted TEs) | genome mode (`-genome`, assembly) |
|---|---|---|
| **Intact / low-divergence sequence** | **pyHMMER** : fastest | **BATH** (only engine) |
| **Degraded / frameshifted copies** | **BATH** : slower than pyHMMER, more sensitive | **BATH** (only engine) : both faster *and* more sensitive |

Genome mode offers no engine choice: BATH is both faster and more sensitive on long sequences, so
the pyHMMER genome path was removed rather than kept as a slower, less sensitive alternative.

### pyHMMER

Single-pass `--nobias` search against all models via pyHMMER, in-process: `hmmsearch` for
amino-acid databases, `nhmmer` for DNA databases (`sine`, `sine-animals`, Dfam).

- Model-cost-aware parallel load balancing chooses between pyHMMER's `queries` and `targets`
  parallelization per model bin, reaching near-full CPU utilization (see
  [Parallel strategy](#parallel-strategy)).
- Results are written to a SQLite database, so filtering and re-analysis do not require re-running
  the search.

### DNA databases (nhmmer)

DNA profile databases are searched with `nhmmer`, not `hmmsearch`. `hmmsearch` only scores the
strand it is handed, so it misses every minus-strand copy, and pyHMMER rejects sequences over
100k residues outright. `nhmmer` scans both strands in one pass and handles long targets.

On 5 Mb of rice against `sine`, nhmmer records hits on both strands (15,992 `+` / 15,821 `-`)
where hmmsearch records no strand at all, and classifies 50 windows against hmmsearch's 37. The
35 windows both engines classify are almost all `+` or unstranded; every one of the 15 windows
only nhmmer recovers is on the minus strand, so the gain is the strand hmmsearch cannot see
rather than a looser threshold.

The comparison above used the earlier `hmmsearch` path, which has since been removed.

**Known limitation — hit filters on the DNA path are inherited from the protein path.** DNA hits
are filtered with the same thresholds as protein domains (coverage ≥ 20%, E-value ≤ 1e-3,
accuracy ≥ 0.5, normalized score ≥ 0.1). Two consequences on the rice/`sine` fixture:

- **Two of the four filters never discriminate.** Accuracy and coverage reject nothing: the same
  50 windows are classified whether the cutoffs are at their defaults or at zero. Only E-value
  and normalized score bind.
- **The E-value cutoff is not comparable between engines.** nhmmer scores against a long-target
  search space spanning both strands, hmmsearch against a per-sequence protein-style one, so the
  same alignment gets very different E-values — for `SHANSINE_MT` on one rice window, 5.5e-05
  under hmmsearch versus 0.0014 under nhmmer. A single `1e-3` cutoff is therefore stricter for
  nhmmer than for hmmsearch, and the two windows hmmsearch classifies that nhmmer does not are
  both boundary misses of this kind, not detection failures — nhmmer finds hundreds of raw hits
  in each.

The practical effect is bounded: dropping the E-value filter altogether raises the count from 50
to 59, so no threshold choice recovers much more. SINE-specific filters are not implemented.

### BATH

The `bath` stage runs `bathsearch --fs` for amino-acid databases. BATH aligns protein pHMMs
directly against raw nucleotide input, allowing the alignment to change frame at indels and read
through stop codons under a frameshift penalty, so a domain interrupted by a frameshift is
recovered as a single hit with its true extent. No six-frame translation is performed.

- HMMER3 databases are converted on first use (`bathconvert`, with `--fs` on BATH builds that support it) and cached, once per run before the search is split across processors. Set `TESORTER2_CACHE` to choose the cache directory.
- Hits are normalized into the same internal schema as the HMMER path, so classification,
  reconciliation and BLAST pass-2 are unchanged.
- BATH tblout reports no per-residue posterior probability, so `acc` is set to `1.0` and the
  minimum-accuracy filter is a no-op for BATH hits.
- Minus-strand coordinates are normalized to ascending order, with strand encoded in the target
  suffix, matching the HMMER convention.
- The only protein engine in genome mode.

### Genome mode (`-genome`)

Treats the input as whole-genome sequence rather than pre-extracted elements: detects TE protein
domains throughout, classifies **each domain individually**, resolves overlapping features, and
emits a domain-level GFF3 plus a summary table. It does not produce a per-element `.cls.tsv` and
does not run BLAST pass-2, matching TEsorter's `-genome` behaviour.

Protein domains are found with **BATH only**. `bathsearch --fs` runs directly on the nucleotide
chromosomes; the tblout already reports nucleotide coordinates and strand, so no six-frame translation
or coordinate back-mapping happens at all. BATH additionally streams long targets in ~0.25 Mb
blocks overlapped by the maximum expected hit length, reconciling boundary duplicates internally.

DNA databases (`sine`) are searched with nhmmer, which streams long targets the same way. Neither
engine needs the genome windowed, so it is searched as given. A DNA-only run is supported and needs
no BATH.

Overlaps are resolved within each feature type, so a BATH protein-domain feature (`CDS`) and an
nhmmer element feature (e.g. `SINE_element`) at the same locus coexist as independent layers rather
than evicting one another.

Requires at least one database.

## Multiple databases

**Most bundled databases are searched by default.** Omitting `-db` searches rexdb, gydb, line, tir,
sine, sine-animals, pfam-te and dfam-core; pass `-db rexdb,gydb` to restrict the set, or
`-db all` to force every alias including the ones that are off by default (the
deeper Dfam collections). Each primary database classifies independently and reconciliation resolves
them afterwards, so more databases means more evidence rather than more ambiguity — at
proportionally more compute.

### Subordinate databases

The four Dfam collections (`dfam-core`, `dfam-extended`, `dfam-deep`, `dfam-complete`) are a
**fallback layer, not competitors**. They run after every primary database has classified what it
can, against only the sequences that were left unclassified, and they are reconciled among
themselves rather than entering the cross-database vote. The protein cascade is the stronger
evidence where it speaks, so a DNA consensus model is never given the chance to outvote it — it only
fills gaps. The layer runs before BLAST pass-2, so the weaker nucleotide-similarity pass sees only
what profile evidence could not reach.

The aliases are **hierarchical**: `-db dfam-deep` searches `dfam-core`, `dfam-extended` and
`dfam-deep`. The files themselves are disjoint, so nothing is searched twice.

Each collection is the previous one plus the next tranche of models, ordered by measured cost per
correctly recovered sequence, so each step costs roughly an order of magnitude more per label than
the one before it. `dfam-core` is on by default; the rest are opt-in.

The three reconciliation stages:

1. **Independent classification**: each database classifies every element on its own, emitting a
   per-database `{prefix}.{db}.cls.tsv` in native TEsorter format.
2. **Name harmonization**: per-database calls are projected onto a unified taxonomy that collapses
   superfamily synonyms (e.g. `Pao` → `Bel-Pao`) and unifies clade names. Lineages with no
   established equivalent are kept distinct to avoid spurious agreement.
3. **Scope-aware hierarchical reconciliation**: a hierarchical vote at Order → Superfamily →
   Clade. At each level candidate labels are weighted by the summed normalized domain score per
   database, and only entries consistent with the winning label advance. The superfamily level is
   *scope-aware*: a database may only elect a superfamily it models with at least 2 clades.
Every per-database call is retained in the `SecondaryHits` column as
`db:order/superfamily/clade=score`, in descending order of evidence — always on the database's
**native** names, so the audit trail survives harmonization. Use `--compat-tesorter-output`
for the original 7-column format.

> **Compatibility contract.** The per-database `{prefix}.{db}.cls.tsv` is **never harmonized** and
> always carries exactly TEsorter's original 7 columns. Downstream pipelines — EDTA reads
> `*.{db}db.cls.tsv` positionally as `(id, Order, Superfamily)` — depend on this, so harmonized
> names and any added columns belong in the combined `{prefix}.cls.tsv` only. Note that
> harmonization rewrites `Order` and `Superfamily` too, not just `Clade` (`Pao`→`Bel-Pao`,
> `pararetrovirus`→`LTR`/`Caulimoviridae`), which is precisely why it must not reach the per-database
> file. If you ever feed the *combined* file to EDTA, extend the `%lib` hash in its
> `cleanup_misclas.pl` first — it knows `pararetrovirus` but not `Caulimoviridae` or `Bel-Pao`.

Harmonization is gated on at least two databases voting on an element: a single-database run keeps
native names throughout, so `-db rexdb` alone is unaffected. The two tables driving stages 2 and 3
are `tesorter2/database/clade_harmonization.tsv` (mapping) and `clade_scope.tsv` (which superfamily
each database can resolve at lineage level). The REXdb↔GyDB lineage equivalences — `Ale`/`Retrofit`,
`Ivana`/`Oryco`, `Tekay`/`Del`, `SIRE`, `Tork`, `Reina`, `CRM`, `Galadriel`, `Athila` — are asserted
from [Neumann et al. 2019](https://doi.org/10.1186/s13100-018-0144-1) (*Mob DNA* 10:1), not fitted
to any dataset. Both tables are read from the resolved `--db-dir` first and fall back to the
packaged copies, so a custom database collection still harmonizes.

Where one database resolves *below* the level both can express, the combined file adds a `Lineage`
column rather than forking the clade name. REXdb splits the Tat group into `Ogre`/`Retand`/`TatI-III`
while GyDB has a single `tat`: both are reported as `Clade=Tat`, and REXdb's finer call lands in
`Lineage`. Lineage is chosen among the databases that actually resolve one, so the finer call
survives even when the coarser database wins the score vote. (`Tatius` is not part of the group —
the REXdb HMM places it at `OTA/Tatius`, a sibling of `Tat`.)

## Sequence Ontology

`Order/Superfamily/Clade` is TEsorter's vocabulary, not a standard one. Every classification is
also resolved to a [Sequence Ontology](http://www.sequenceontology.org) term, so results are
comparable with other annotation tools:

- `{prefix}.cls.tsv` gains `SO_name` and `SO_ID` columns (e.g. `Copia_LTR_retrotransposon`,
  `SO:0002264`).
- Genome-mode GFF3 features carry `Ontology_term=SO:...` plus `so_name=`. The feature type stays
  `CDS`: these features are protein domains, not elements, so typing one as
  `Gypsy_LTR_retrotransposon` would assert the domain *is* the retrotransposon.

The mapping authority is [EDTA's `TE_Sequence_Ontology.txt`](https://github.com/oushujun/EDTA/blob/master/bin/TE_Sequence_Ontology.txt),
bundled in `tesorter2/data/`. All 59 Order/Superfamily labels the bundled databases can emit
resolve to a specific SO term; nothing falls back to the generic `repeat_region`.

For lineages with no SO term of their own, EDTA files a descriptive name under a generic
accession (`CR1_LINE_retrotransposon` is not a real SO term; its `SO:0000194` is
`LINE_element`'s). `SO_name` keeps EDTA's name for interoperability, while anything written as an
ontology term resolves to the real one.

`--compat-tesorter-output` suppresses the SO columns, keeping the original 7-column format.

## Clade voting

Within each database the winning clade is chosen by a **score-weighted vote**: each domain
contributes a length-normalized score (`dom_score / model_len`) to its clade, and the highest
summed score wins:

$$\hat{c} = \arg\max_{c} \sum_{d \in D_c} \frac{\text{score}_d}{\text{len}_d}$$

Normalizing by model length controls for the tendency of longer profiles to accumulate higher raw
scores, making domains comparable. Ties fall through to the existing mixture/completeness rules;
Order and Superfamily are inherited from the winning clade.

This replaces TEsorter's raw domain-count plurality, which breaks ties by position in an internal
collection rather than by any biological signal, causing sibling-clade swaps (Reina↔Tekay,
Ale↔Alesia) to flip depending on search engine. Pass `--compat-tesorter-voting` to restore the
original behaviour.

## CLI reference

```
TEsorter2 <fasta> [options]               # also installed as tesorter2
TEsorter  <fasta> [options]               # TEsorter's own pipeline
```

The command line is TEsorter's, plus TEsorter2's own options. No mode has to be named: element
classification is the default and `-genome` switches to genome mode. Output files go beside
`-pre`, which defaults to `{input}.{db}` with one database (TEsorter's naming) and to
`{input}.TEsorter2` with several, in the working directory.

TEsorter's options:

| Flag | Default | Description |
|---|---|---|
| `-db`, `--hmm-database` | all bundled except the deeper Dfam (`TEsorter`: `rexdb`) | Database name, HMM file, comma-separated list, or `all`. TEsorter's names are accepted: `rexdb-plant`, `rexdb-metazoa` and the v3 names search the bundled REXdb (plant v4.0 + metazoa v3.1), `rexdb-pnas` is `tir`, `rexdb-line` is `line` |
| `--db-hmm` | none | Custom HMM file, used instead of `-db` |
| `-st`, `--seq-type` | `nucl` | `prot` for amino-acid input: no six-frame translation, no BATH, DNA databases skipped |
| `-pre`, `--prefix` | `{input}.{db}` / `{input}.TEsorter2` | Output prefix, may include a directory |
| `-p`, `--processors` | `4` | Processors |
| `-tmp`, `--tmp-dir` | `./tmp-{random id}` | Temporary files; removed at the end unless `-nocln` |
| `-cov`, `-eval`, `-prob`, `-score` | `20`, `1e-3`, `0.5`, `0.1` | Domain filters: HMM coverage (%), E-value, posterior probability, score per HMM position |
| `-mask soft hard` | none | Write `{prefix}.softmasked` / `.hardmasked` with the TE domains masked. Needs domain coordinates: `-genome`, or the `TEsorter` command |
| `-nocln`, `--no-cleanup` | cleanup **on** | Keep intermediates (partitions, translations, cascade handoffs, pass-2 scratch). A failed run always keeps them |
| `-dp2`, `--disable-pass2` | off | Skip the pass-2 similarity search (`TEsorter` skips it for `-st prot`, as TEsorter did) |
| `-rule`, `--pass2-rule` | `80-80-80` | Pass-2 threshold, identity-coverage-length |
| `-nolib`, `-norc` | off | No RepeatMasker library; do not reverse-complement minus-strand sequences in it |
| `-genome` | off | Genome mode |
| `-fw`, `-win_size`, `-win_ovl` | | Accepted and ignored: there is no hmmscan output to reuse, and genome mode searches each chromosome whole |
| `-cite`, `-v` | | Print the citation / version and exit |

TEsorter2's options:

| Flag | Default | Description |
|---|---|---|
| `-o`, `--outdir` | directory of `-pre` | Output directory; a relative `-pre` is placed inside it |
| `--db-dir` | bundled | Directory holding the HMM databases (or `TESORTER2_DB`) |
| `--mask-stops` | off | Translate stop codons as `X` instead of `*`, letting profiles align through premature stops |
| `--stages` | `nail,hmmer,bath` | Engine cascade for amino-acid databases (see above) |
| `--no-partition` | partitioning **on** | Search the input as one unit instead of one balanced group per processor (partitioning is 4.6x faster on 135k sequences over 64 cores) |
| `--memory-db` | off | Build the results database in RAM and copy it out once at the end |
| `--extend DB` | off | Add the `-db` databases to a finished run's results database; vote, pass-2 and outputs are redone |
| `--emit-timing` | off | Write `cascade_timing.tsv` and `stage_timing.tsv`. Diagnostic |
| `--pass2-aligner` | `blast` | `blast` (blastn, or blastp for `-st prot`) or `minimap2` |
| `--blast-task` | `megablast` | `dc-megablast` is slower but more sensitive to diverged matches |
| `--pass2-classified-fasta` | none | Extra classified elements (`>id#Order/Superfamily/Clade`) for the pass-2 target pool |
| `--minimap2-extra` | empty | Extra flags passed to minimap2 |
| `--min-clade-delta BITS` | `0` | Demote a clade call to `mixture` when a score offset below BITS could overturn it |
| `--emit-clade-delta` | off | Add a `CladeDelta` column to the per-database `.cls.tsv`. Diagnostic |
| `--compat-tesorter-voting` | off | Raw domain-count clade vote (TEsorter behaviour) |
| `--compat-tesorter-rounding` | off | Round normalized scores to 2 dp before filtering (replicates a TEsorter bug) |
| `--compat-tesorter-output` | off | Combined `.cls.tsv` in TEsorter's 7-column format |
| `--no-tesorter-outputs` | off | Skip the per-database RepeatMasker libraries |

The 2.0.0 spellings still work for existing scripts but are not listed in `-h`: the subcommand
words (`tesorter2 sequences lib.fa`, `tesorter2 genome asm.fa`), `-d`/`--database` for `-db`,
`--max-search` for `-db all`, and `--genome` for `-genome`.

Genome mode uses `-db`, `--db-hmm`, `-pre`, `-o`, `-tmp`, `--db-dir`, `-p`, the
domain filters and `-mask`; the element-mode options are ignored there.

## Pass-2 similarity search

Sequences that no database classified are compared with the classified ones. With the default
`blast` aligner, a hit passing the I-C-L rule (identity, query coverage, alignment length) assigns
Order and Superfamily with `clade=unknown`. With `--pass2-aligner minimap2`, pass-2 runs

```
minimap2 -x asm20 --rmq=no --no-long-join \
    -k 10 -w 10 -r 500,20000 -g 500 -p 0.3 -N 100 -m 30 \
    -t NCPU -K 1G --seed 11 --paf-no-hit \
    -o pass2.paf  TARGET.fa  QUERY.fa
```

and the PAF is collapsed by `classify_ltr_paf_fast` to one row per query; a query inherits the best
target's full classification iff pid ≥ I/100, eff_qcov ≥ C/100 and eff_tcov ≥ C/100 (L is not
used). Benchmark at `70-70-70`: F1 ≈ 0.895, accuracy ≈ 0.943, precision ≈ 0.866, recall ≈ 0.926,
MCC ≈ 0.857. The `blast_hits` table of the results database stores one row per query (best target):
`qseqid, sseqid, pident, qcovs, tcovs, passes_rule, classified_by`.

## TEsorter compatibility

`TEsorter2` takes TEsorter's command line, so a TEsorter command keeps working when only the
executable name changes; it runs the TEsorter2 search with TEsorter's file names. For workflows
built around TEsorter 1, `TEsorter` runs TEsorter's own pipeline: one database (rexdb by default), HMMER on the
six-frame translation, TEsorter's filters and count-based clade voting:

```bash
TEsorter input.fasta -db rexdb -p 16 -pre out
```

Every TEsorter 1 option is accepted (table above). Output file names follow TEsorter
(`{input}.{db}.cls.tsv`, `.cls.lib`, `.cls.pep`, `.dom.gff3`, `.dom.tsv`, `.dom.faa`, and
`.softmasked`/`.hardmasked` with `-mask`), so EDTA can call it unchanged. With `-genome` it runs
TEsorter2's genome mode on the one database, under TEsorter's genome-mode names.

## Architecture

### Core

- **`cli.py`** — the command line shared by `TEsorter2` and `TEsorter`
- **`pipeline.py`** — search orchestration
- **`tesorter_compat.py`** — TEsorter's own pipeline (the `TEsorter` command)
- **`search.py`** — HMM search engine with balanced parallelism
- **`sequence.py`** — FASTA ingestion (pyfastx) and six-frame translation (pyhmmer)
- **`hmm.py`** — HMM loading, alphabet detection, optimized profile construction
- **`bath_search.py`** — BATH engine: conversion, `bathsearch --fs` invocation, hit normalization
- **`nail_search.py`**, **`nhmmer_search.py`**, **`hierarchical_search.py`** — nail and nhmmer engines; the staged cascade
- **`partition.py`** — splitting the input across processors and merging the results
- **`genome.py`** — Genome mode: whole-chromosome search, per-domain classification, overlap resolution, GFF3/summary
- **`results.py`** — SQLite persistence with pre-parsed columns (base_seq, strand, frame, domain_type)
- **`deconflict.py`** — numpy-based hit deconfliction and parameterized filtering
### Classification and post-processing

- **`classifier.py`** — config-driven classification from domain hits; per-database domain
  remapping, overlap-aware deconfliction, order/superfamily/clade assignment
- **`blast_pass2.py`**, **`minimap.py`**, **`classify_ltr_paf_fast.py`** — pass-2 with BLAST or minimap2, cross-database target pooling
- **`extend.py`** — `--extend`: adding databases to a finished run
- **`tesorter_compat.py`** — the `TEsorter` command line
- **`tesorter_output.py`**, **`emit.py`**, **`id_registry.py`** — output formatting and identifier bookkeeping

## Parallel strategy

**Default search.** pyHMMER exposes two C-level parallelization schemes: `queries`, where each
thread takes one HMM and searches it against all sequences, and `targets`, where one sequence is
searched against models in parallel. `queries` is inherently more efficient unless there are few
models.

HMM runtime scales roughly with *M²*, so a single long model can dominate: under `queries`
parallelism every other model finishes quickly while the long one runs single-threaded.

TEsorter2 precomputes each model's expected cost (*M²*) and bins them: **small** models
(cost ≤ 75th percentile + 2×IQR for that database) run in `queries` mode; **large** models run in
`targets` mode. Sequences are reused from the same in-memory object across both searches, so the
split is essentially free, yielding near-full CPU utilization in the most efficient mode available
for each model class.
