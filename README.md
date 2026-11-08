# TEsorter2

Fast, divergence-robust classification of transposable elements (TEs) from their protein domains.

TEsorter2 is the successor of [TEsorter](https://github.com/zhangrengang/TEsorter). It keeps
TEsorter's classification scheme and its command line, and adds:

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
- **Standard vocabulary.** Every call also carries its [Sequence Ontology](http://www.sequenceontology.org) term.
- **Drop-in replacement.** TEsorter2 takes TEsorter's command line unchanged (`-db`, `-pre`, `-p`,
  `-cov`, `-genome`, ...): replacing `TEsorter` with `TEsorter2` in a command is enough. The
  `TEsorter` command is also installed and runs TEsorter's own pipeline.

## Installation

### conda

```bash
conda install -c conda-forge -c bioconda tesorter
```

This installs TEsorter2 with every search engine (HMMER, nail with MMseqs2, BATH) plus BLAST+
and minimap2. The HMM databases are bundled.

### Container

A Docker / Apptainer image with every engine is described in
[`container/README.md`](container/README.md):

```bash
docker build -f container/Dockerfile -t tesorter2:2.0.0 .
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data tesorter2:2.0.0 TEsorter2 TEs.fa -p 8
```

### From source

```bash
git clone https://github.com/KGerhardt/TEsorter2.git
cd TEsorter2
conda env create -f environment.yml
conda activate tesorter2
```

or, with HMMER and BLAST+ already on `PATH`, `pip install git+https://github.com/KGerhardt/TEsorter2.git`.

### BATH outside conda

With a pip or source install, BATH has to be added by hand (it is optional when classifying TE
sequences, where the BATH stage is skipped with a warning, and required for `-genome`). Download a prebuilt release from the [BATH releases page](https://github.com/TravisWheelerLab/BATH/releases)
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

The command line is TEsorter's: an input FASTA, plus options. There is no mode to choose;
`-genome` switches to whole-genome annotation, as in TEsorter.

### Classify a TE library

```bash
TEsorter2 TElib.fa -p 16
```

With no `-db`, every default database is searched (rexdb, gydb, line, tir, sine, sine-animals,
pfam-te, dfam-core). Output files go to the working directory, named `TElib.fa.TEsorter2.*`
(change with `-pre`). The main file is `TElib.fa.TEsorter2.cls.tsv`, one row per classified sequence
(tab-separated; the last columns are shortened here):

```
#TE                    Order  Superfamily  Clade      Complete  Strand  Domains                SecondaryHits                                      SO_name                    SO_ID       Lineage
Os1468_INT#LTR/Gypsy   LTR    Gypsy        Tekay/Del  yes       +       GAG|del AP|del RT|del…  gydb:LTR/Gypsy/del=8.984;rexdb:LTR/Gypsy/Tekay=8.661  Gypsy_LTR_retrotransposon  SO:0002265  .
```

`SecondaryHits` lists every database's own call with its score, strongest first. Each database
also gets its own `TElib.fa.TEsorter2.{db}.cls.tsv` in TEsorter's original format, and
`TElib.fa.TEsorter2.cls.lib` is the RepeatMasker library.

### One database, as in TEsorter

```bash
TEsorter2 TElib.fa -db rexdb-plant -p 16
```

This writes `TElib.fa.rexdb-plant.cls.tsv` and `TElib.fa.rexdb-plant.cls.lib`, the names TEsorter
uses, with the calls of the TEsorter2 engine cascade. TEsorter's database names (`rexdb`,
`rexdb-plant`, `rexdb-metazoa`, `gydb`, `rexdb-pnas`, `rexdb-line`, `sine`) work, as do
TEsorter2's; `-db` also takes a comma-separated list, an HMM file, or `all`.

### Annotate TE domains in a genome

```bash
TEsorter2 assembly.fa -genome -p 16
```

This writes `assembly.fa.TEsorter2.dom.gff3` (one feature per classified domain, with its
classification and SO term) and `assembly.fa.TEsorter2.genome.summary.tsv` (counts per Order /
Superfamily / Clade). Genome mode needs BATH for the protein databases. `-mask soft|hard` also
writes the masked assembly.

### TEsorter's own pipeline

```bash
TEsorter TElib.fa -db rexdb-plant -p 16
```

TEsorter's search (one database, rexdb by default, HMMER only), filters and count-based clade vote,
with every TEsorter file (`.cls.tsv`, `.cls.lib`, `.cls.pep`, `.dom.gff3`, `.dom.tsv`, `.dom.faa`).
Use it where downstream tools expect TEsorter 1 results exactly, as EDTA does.

## Common options

```bash
# Choose databases
TEsorter2 TElib.fa -db rexdb,gydb

# Add the deeper Dfam collections (dfam-extended < dfam-deep < dfam-complete), or everything
TEsorter2 TElib.fa -db rexdb,dfam-deep
TEsorter2 TElib.fa -db all

# Choose engines: fastest (HMMER only), or BATH only, or the full default cascade
TEsorter2 TElib.fa --stages hmmer
TEsorter2 TElib.fa --stages bath
TEsorter2 TElib.fa --stages nail,hmmer,bath

# Protein input
TEsorter2 proteins.fa -st prot -db rexdb

# TEsorter's domain filters (defaults shown)
TEsorter2 TElib.fa -cov 20 -eval 1e-3 -prob 0.5 -score 0.1

# Pass-2 (similarity rescue of unclassified sequences): skip it, or use minimap2
TEsorter2 TElib.fa -dp2
TEsorter2 TElib.fa --pass2-aligner minimap2 -rule 80-80-80

# Add a database to a finished run without searching everything again
TEsorter2 TElib.fa --extend TElib.fa.TEsorter2.db -db dfam-extended

# Output name and location; temporary files
TEsorter2 TElib.fa -pre results/mylib -tmp /scratch/tmp
```

Run `TEsorter2 -h` for the full list.

## Choosing an engine

| Your sequences | Suggested run |
|---|---|
| Recent, intact elements (e.g. a fresh EDTA / RepeatModeler library) | default, or `--stages hmmer` for maximum speed |
| Old, diverged or fragmented copies | default cascade with BATH installed |
| Whole assembly | `-genome` (BATH) |

In the default cascade each engine only searches what the previous ones left unclassified, so
BATH, the slowest and most sensitive, sees a small remainder. DNA databases (SINE, Dfam) are
always searched with nhmmer.

## Output files

`{prefix}` is `-pre`, by default `{input}.{db}` with one database and `{input}.TEsorter2` with several.

| File | Content |
|---|---|
| `{prefix}.cls.tsv` | Combined call per sequence across databases and pass-2, with `SecondaryHits`, `SO_name`, `SO_ID` |
| `{prefix}.cls.lib` | RepeatMasker library of the combined calls (`>id#Order/Superfamily/Clade`) |
| `{prefix}.{db}.cls.tsv` | Several databases: per-database call in TEsorter's original 7-column format (never harmonized) |
| `{prefix}.{db}.cls.lib` | Several databases: per-database RepeatMasker library |
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
