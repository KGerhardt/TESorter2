# Bundled database licenses

TEsorter2 ships several HMM databases inside the package (`tesorter2/database/`).
These databases were **not** produced by the TEsorter2 authors; each carries its
own upstream license and citation requirement, listed below. The TEsorter2 code
is GPL-3.0-or-later (see the top-level `LICENSE`), but that license does **not**
extend to the third-party data in this directory — the terms below govern the
`.hmm` files and their metadata.

If you use TEsorter2 in published work, cite both TEsorter2 and the databases you
actually ran against (the alias you pass to `-d` selects the database).

| File(s) | Database | Alias | License | Cite |
|---|---|---|---|---|
| `REXdb_protein_database_viridiplantae_v4.0_plus_metazoa_v3.1.hmm` | REXdb | `rexdb` | CC BY 4.0 | Neumann et al. 2019 |
| `GyDB2.hmm`, `GyDB2.hmm.info` | GyDB 2.0 | `gydb` | CC BY (Attribution) | Llorens et al. 2011 — **modified, see notice** |
| `AnnoSINE_core.hmm` | AnnoSINE | `sine` | MIT | Li et al. 2022 |
| `Kapitonov_et_al.GENE.LINE.hmm` | Kapitonov LINE | `line` | see note | Kapitonov et al. 2009 |
| `Yuan_and_Wessler.PNAS.TIR.hmm` | Yuan & Wessler TIR | `tir` | see note | Yuan & Wessler 2011 |


---

## REXdb — CC BY 4.0

Protein-domain profiles for LTR-retrotransposons (Viridiplantae v4.0 + Metazoa v3.1).

- **License:** Creative Commons Attribution 4.0 International (CC BY 4.0), as
  published on the RepeatExplorer Zenodo record.
- **Source:** <https://github.com/repeatexplorer/rexdb> · <http://repeatexplorer.org/?page_id=918>
- **Cite:** Neumann P, Novák P, Hoštáková N, Macas J (2019). *Systematic survey of
  plant LTR-retrotransposons elucidates phylogenetic relationships of their
  polyprotein domains and provides a reference for element classification.*
  Mobile DNA 10:1. doi:10.1186/s13100-018-0144-1

## GyDB 2.0 — Creative Commons Attribution

Profiles for the Gypsy Database of mobile genetic elements (release 2.0).

- **License:** Creative Commons Attribution License, as stated on the GyDB site
  ("distributed under the terms of the Creative Commons Attribution License,
  which allows unrestricted use, distribution, and reproduction in any medium,
  provided the original work is properly cited").
- **Source:** <https://gydb.org/index.php/The_Gypsy_Database:About>
- **Cite:** Llorens C, Futami R, Covelli L, et al. (2011). *The Gypsy Database
  (GyDB) of mobile genetic elements: release 2.0.* Nucleic Acids Research
  39(suppl_1):D70–D74. doi:10.1093/nar/gkq1061

### Modification notice (CC BY requires changes to be indicated)

`GyDB2.hmm` as shipped here **has been modified** from the file distributed by
GyDB. The change adds the `COMPO` record — the model's average residue
composition — to all 314 profiles, and updates the format header from
`HMMER3/f [3.1b2 | February 2015]` to `HMMER3/f [3.4 | Aug 2023]`.

**Why.** GyDB's profiles were built in 2009 with HMMER 3.1b2, before `COMPO`
became a standard part of the HMMER3 save format. Every other database bundled
here carries it (REXdb 266/266, LINE 28/28, TIR 17/17); GyDB2 had 0/314. HMMER
and pyHMMER tolerate the omission by computing composition internally, so it is
invisible to them, but tools that require a well-formed HMMER3 record reject the
file outright — `nail` fails with *"failed to write mmseqs query DB / missing:
match emissions 1"* and cannot search GyDB at all.

**What was not changed.** No emission or transition probability was altered, and
no model was rebuilt. The composition is *derived* from each model's existing
match emissions by HMMER's own `p7_hmm_SetComposition()` (via pyHMMER's
`HMM.set_composition()`), not re-estimated from sequence data. Rebuilding from
the original GyDB alignments with a current `hmmbuild` was rejected deliberately:
priors and effective-sequence weighting changed between 3.1b2 and 3.4, so that
route would produce different models and shift existing results. All 314 model
names and lengths are unchanged.

**Verified equivalent.** On 500 LTR RepBase elements (25,470 hits), the original
and modified files give byte-identical hit lists — target, profile, domain score,
HMM bounds, i-Evalue — and identical classifications.

**Reproduce:**

```python
import pyhmmer.plan7 as p7
with p7.HMMFile(original) as fh, open(modified, "wb") as out:
    for hmm in fh:
        if hmm.composition is None:
            hmm.set_composition()
        hmm.write(out)
```

## AnnoSINE — MIT

SINE profiles (`AnnoSINE_core.hmm`), resolved by the `sine` alias.

- **License:** MIT (per the AnnoSINE repository).
- **Source:** <https://github.com/yangli557/AnnoSINE>
- **Cite:** Li Y, Jiang N, Sun Y (2022). *AnnoSINE: a short interspersed nuclear
  elements annotation tool for plant genomes.* Plant Physiology 188(2):955–970.
  doi:10.1093/plphys/kiab524

## Kapitonov LINE — no separate license; redistributed via TEsorter (GPL-3.0)

HMM profiles of non-LTR retrotransposon (LINE) RT/domain clades. The profiles are
distributed as part of TEsorter (GPL-3.0), which is the immediate upstream of
these files; no standalone license is published for the profiles themselves.

- **Upstream package:** TEsorter — <https://github.com/zhangrengang/TEsorter> (GPL-3.0)
- **Cite:** Kapitonov VV, Tempel S, Jurka J (2009). *Simple and fast
  classification of non-LTR retrotransposons based on phylogeny of their RT
  domain protein sequences.* Gene 448(2):207–213. doi:10.1016/j.gene.2009.07.019

## Yuan & Wessler TIR — no separate license; redistributed via TEsorter (GPL-3.0)

HMM profiles of the catalytic domains of cut-and-paste (TIR) transposase
superfamilies. As with the LINE profiles, these are redistributed from TEsorter
(GPL-3.0); no standalone license is published for the profiles themselves.

- **Upstream package:** TEsorter — <https://github.com/zhangrengang/TEsorter> (GPL-3.0)
- **Cite:** Yuan YW, Wessler SR (2011). *The catalytic domain of all eukaryotic
  cut-and-paste transposase superfamilies.* PNAS 108(19):7884–7889.
  doi:10.1073/pnas.1104208108

---

*Verified 2026-07-17. If an upstream project changes its license, the terms in
effect when you obtained the database apply; re-check the source links above.*
