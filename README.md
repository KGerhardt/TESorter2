# TEsorter2

Fast, divergence-robust classification of transposable elements (TEs) from their protein domains.

TEsorter2 is the successor of [TEsorter](https://github.com/zhangrengang/TEsorter). It keeps
TEsorter's classification scheme (Order / Superfamily / Clade from conserved TE protein domains)
and its command line, and adds:

- **Speed.** Searches run in-process through [pyHMMER](https://pyhmmer.readthedocs.io/) with
  cost-aware load balancing, and the input is split across processors. On 108k RepBase elements
  TEsorter2 is **5.4× faster** than TEsorter; on a whole rice genome, **4.9× faster**.
- **Sensitivity on old, degraded copies.** A staged cascade (`nail → hmmer → BATH`) sends what the
  fast engines cannot classify to [BATH](https://github.com/TravisWheelerLab/BATH), which aligns
  protein profiles directly to DNA and tolerates frameshifts and stop codons. Domains broken by
  indels are recovered as a single hit.
- **Many databases in one run, reconciled.** REXdb, GyDB, LINE, TIR, SINE, a TE-focused Pfam set
  and Dfam are bundled. Each classifies independently and a score-weighted, scope-aware vote
  reconciles them, keeping every per-database call for auditing.
- **Reproducible calls.** The clade is chosen by a score-weighted vote instead of a raw domain count
  whose ties depended on internal ordering.
- **Genome mode.** Annotates TE domains along a whole assembly and writes a GFF3.
- **Standard vocabulary.** Every call also carries its [Sequence Ontology](http://www.sequenceontology.org)
  term.
- **Drop-in replacement.** The `TEsorter` command accepts TEsorter's arguments and writes the same
  files, so pipelines such as [EDTA](https://github.com/oushujun/EDTA) keep working unchanged.

## Installation

### conda

```bash
conda install -c conda-forge -c bioconda tesorter
```

This installs TEsorter2 with HMMER, BLAST+, minimap2, nail and MMseqs2, and provides both the
`tesorter2` and the `TEsorter` commands. The HMM databases are bundled: there is nothing to
download. BATH is not on conda yet; see below.

### Container (includes BATH)

A Docker / Apptainer image with every engine, BATH included, is described in
[`container/README.md`](container/README.md):

```bash
docker build -f container/Dockerfile -t tesorter2:2.0.0 .
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data tesorter2:2.0.0 tesorter2 sequences TEs.fa -p 8
```

### From source

```bash
git clone https://github.com/KGerhardt/TEsorter2.git
cd TEsorter2
conda env create -f environment.yml
conda activate tesorter2
```

or, with HMMER and BLAST+ already on `PATH`, `pip install git+https://github.com/KGerhardt/TEsorter2.git`.

### Adding BATH (optional for `sequences`, required for `genome`)

Without BATH, `tesorter2 sequences` simply skips the BATH stage and says so in the log. To enable
it, download a prebuilt release from the [BATH releases page](https://github.com/TravisWheelerLab/BATH/releases)
or build it from source:

```bash
git clone https://github.com/TravisWheelerLab/BATH.git
cd BATH
git clone -b BATH https://github.com/TravisWheelerLab/easel
autoconf && ./configure && make
export PATH=$PWD/src:$PATH        # or: export BATH_BIN_DIR=$PWD/src
```

TEsorter2 converts each database for BATH on first use and caches it.

## Quick start

### Classify a TE library

```bash
tesorter2 sequences TElib.fa -p 16
```

Results go to `TElib.fa.TEsorter2/` (change with `-o`). The main file is `TElib.fa.cls.tsv`, one
row per classified sequence (tab-separated; the last columns are shortened here):

```
#TE                    Order  Superfamily  Clade      Complete  Strand  Domains                SecondaryHits                                      SO_name                    SO_ID       Lineage
Os1468_INT#LTR/Gypsy   LTR    Gypsy        Tekay/Del  yes       +       GAG|del AP|del RT|del…  gydb:LTR/Gypsy/del=8.984;rexdb:LTR/Gypsy/Tekay=8.661  Gypsy_LTR_retrotransposon  SO:0002265  .
```

`SecondaryHits` lists every database's own call with its score, strongest first. Each database
also gets its own `TElib.fa.{db}.cls.tsv` in TEsorter's original format.

### Annotate TE domains in a genome

```bash
tesorter2 genome assembly.fa -p 16
```

This writes, in `assembly.fa.TEsorter2/`, `assembly.fa.dom.gff3` (one feature per classified domain, with its classification and
SO term) and `assembly.fa.genome.summary.tsv` (counts per Order / Superfamily / Clade). Genome mode
needs BATH for the protein databases.

### Use it in place of TEsorter

```bash
TEsorter TElib.fa -db rexdb-plant -p 16
```

Same arguments, same output names (`TElib.fa.rexdb-plant.cls.tsv`, `.cls.lib`, `.dom.gff3`, ...)
as TEsorter 1, with TEsorter's filters and count-based voting. On a rice TE library it agrees with
TEsorter 1 on 99.7% of calls (100% at Order). Existing scripts and EDTA need no changes.

## Common options

```bash
# Choose databases (default: rexdb, gydb, line, tir, sine, sine-animals, pfam-te, dfam-core)
tesorter2 sequences TElib.fa -d rexdb,gydb

# Add the deeper Dfam collections (dfam-extended < dfam-deep < dfam-complete)
tesorter2 sequences TElib.fa -d rexdb,dfam-deep

# Choose engines: fastest (HMMER only), or BATH only, or the full default cascade
tesorter2 sequences TElib.fa --stages hmmer
tesorter2 sequences TElib.fa --stages bath
tesorter2 sequences TElib.fa --stages nail,hmmer,bath

# Protein input
tesorter2 sequences proteins.fa -st prot -d rexdb

# Pass-2 (similarity rescue of unclassified sequences): skip it, or use minimap2
tesorter2 sequences TElib.fa -dp2
tesorter2 sequences TElib.fa --pass2-aligner minimap2 -rule 80-80-80

# Add a database to a finished run without searching everything again
tesorter2 sequences TElib.fa --extend TElib.fa.TEsorter2/TElib.fa.db -d dfam-extended

# Output location and name
tesorter2 sequences TElib.fa -o results/ --prefix mylib
```

Run `tesorter2 sequences -h` or `tesorter2 genome -h` for the full list.

## Choosing an engine

| Your sequences | Suggested run |
|---|---|
| Recent, intact elements (e.g. a fresh EDTA / RepeatModeler library) | default, or `--stages hmmer` for maximum speed |
| Old, diverged or fragmented copies | default cascade with BATH installed |
| Whole assembly | `tesorter2 genome` (BATH) |

In the default cascade each engine only searches what the previous ones left unclassified, so
BATH, the slowest and most sensitive, sees a small remainder. DNA databases (SINE, Dfam) are
always searched with nhmmer.

## Output files

| File | Content |
|---|---|
| `{prefix}.cls.tsv` | Combined call per sequence across databases and pass-2, with `SecondaryHits`, `SO_name`, `SO_ID` |
| `{prefix}.{db}.cls.tsv` | Per-database call in TEsorter's original 7-column format (never harmonized) |
| `{prefix}.{db}.cls.lib` | Per-database RepeatMasker library (`>id#Order/Superfamily`) |
| `{prefix}.db` | SQLite database with every hit and call, for re-analysis without re-searching |
| `{prefix}.dom.gff3` | Genome mode: classified TE domain features |
| `{prefix}.dom.fna` | Genome mode: sequence of each feature |
| `{prefix}.genome.summary.tsv` | Genome mode: Order / Superfamily / Clade tallies |

## Benchmarks

ANVIL CPU node (AMD EPYC 7763), 16 threads, mean of 3 replicates.

| Task | TEsorter | TEsorter2 (pyHMMER) | TEsorter2 (BATH) |
|---|---|---|---|
| 108,318 RepBase LTR/TIR/LINE elements | 2,932 s | **543 s (5.4×)** | 983 s (3.0×) |
| Rice genome, ~375 Mb, REXdb | 2,127 s | — | **431 s (4.9×)** |

## More detail

How the cascade, reconciliation, clade vote, Sequence Ontology mapping, DNA databases and parallel
strategy work, plus the full CLI reference: [`docs/DETAILS.md`](docs/DETAILS.md).

## Citation

If you use TEsorter2, please cite TEsorter2 (manuscript in preparation) and the original TEsorter:

> Zhang R-G, Li G-Y, Wang X-L, *et al.* (2022) TEsorter: an accurate and fast method to classify
> LTR-retrotransposons in plant genomes. *Horticulture Research* 9: uhac017.
> https://doi.org/10.1093/hr/uhac017

Please also cite the databases you search; see
[`tesorter2/database/LICENSES.md`](tesorter2/database/LICENSES.md).

## License

GPL-3.0-or-later. The bundled HMM databases are third-party data with their own licenses
(CC BY 4.0 for REXdb, Creative Commons Attribution for GyDB2, MIT for AnnoSINE; Kapitonov LINE and
Yuan & Wessler TIR are redistributed as in TEsorter). TEsorter2's license does not extend to them;
details in [`tesorter2/database/LICENSES.md`](tesorter2/database/LICENSES.md).
