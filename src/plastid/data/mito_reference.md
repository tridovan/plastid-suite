# The mito reference set — where it came from and what it cannot do

**Why a file exists at all.** A plant mitogenome is annotated by transferring an annotation from
curated relatives (see `docs/tools/pmga.md` for why the alternative — a container driving MAKER —
is not a route this project can take). This is the set that transfer uses. It is an **artifact**:
`mito_reference_proteins.faa` and `mito_reference_rrna.fa` are committed beside this file, and both
are rebuilt from the pinned records by `core/mito_reference.py`.

## Sources — pinned 2026-10-03, NCBI RefSeq, whole mitogenomes

| Accession | Organism | Family | Length (bp) | CDS in the record | sha256 |
|---|---|---|---|---|---|
| `NC_023209.1` | *Salvia miltiorrhiza* | Lamiaceae | 499,236 | 138 | `cea17b5414ce0a61…` |
| `NC_023103.1` | *Ajuga reptans* | Lamiaceae | 352,069 | 30 | `e26e43b35135f2e9…` |
| `NC_049064.1` | *Rotheca serrata* | Lamiaceae | 482,114 | 34 | `27dd00c0f9f5fd4c…` |
| `NC_006581.1` | *Nicotiana tabacum* | Solanaceae | 430,597 | 156 | `ec424b4de5b20e07…` |
| `NC_037304.1` | *Arabidopsis thaliana* | Brassicaceae | 367,808 | 33 | `5b42e0ebb980a282…` |
| `NC_011033.1` | *Oryza sativa* Japonica | Poaceae | 490,520 | 53 | `2e1065c161d421af…` |
| `NC_066984.1` | *Psychotria viridis* | Rubiaceae | 615,370 | 49 | `1db770c1ca267be1…` |
| `NC_057654.1` | *Scyphiphora hydrophyllacea* | Rubiaceae | 354,155 | 37 | `d504cc8b109fddd1…` |

Chosen by the same rule the plastid Reference Panel uses — **same family first**. The two targets
are *Oldenlandia corymbosa* (Rubiaceae) and *Coleus hadiensis* (Lamiaceae), so each has three
relatives of its own family in the set, plus the model species as an outgroup that keeps the core
genes honest. Full sha256 and md5 in `.panel/mito-reference/PINNED.json`; the records themselves
are **not** committed (5.2 MB), and the one-line rebuild is below.

## What the build refuses to include — and says so

Measured 2026-10-03 over the eight records. Every exclusion is returned by the builder and
printed, never dropped in silence:

| Excluded | Count | Why |
|---|---|---|
| anonymous open reading frames (`orf101a`, `orf297a`, …) | **236** | *Nicotiana* annotates 121 of them as CDS, *Salvia* 97, *Oryza* 22. They are not genes with names, and one could outscore the real gene beside it |
| features named by a product sentence (`hypothetical protein`, `LAGLIDADG endonuclease`) | 7 | a sentence is not a symbol |
| rRNA outside its length band | 2 | *Psychotria* calls a **722 bp** feature "5S" and an **1112 bp** one "18S"; plant mito rRNA measures 5S 107–122 · 18S 1692–1935 · 26S 2907–3507 |

**The proteins are the record's own `/translation`, not our translation.** Measured: all **530**
CDS publish one, and **none** declares a `/transl_table`. Self-translating does not merely lose
accuracy — plant mito CDS begin with codons like `ACG` that no table accepts as an initiator, so a
CDS-mode translation raises before producing anything. The curator's protein also gets
trans-splicing and RNA-edited ends right, which this module has no business re-deriving.

## The measured gap — the vocabulary is narrower than the gene set

**31/31 vocabulary loci are covered** by at least one reference protein — that is the artifact's
Done-when, and it holds. But the references carry **23 more genes this toolkit cannot name**:

* **standard plant-mito genes the vocabulary lacks** — `rpl2` · `rpl10` · `rpl16` · `rpl23` ·
  `rps2` · `rps7` · `rps12` · `rps14` · `rps19` · `sdh3`. These are real, conserved, and present
  in most of the eight records. `data/mito_gene_names.yaml` holds 31 loci because it was built
  from one record's own product strings; this is the cross-check that file says it still needs.
* **plastid-derived and horizontally transferred genes** in the mito — `rbcL` · `petG` · `petL` ·
  `atpB` · `atpE` (all in *Salvia*'s record) · `citB` · `epd` · `ltrA` · `rnhA` · `tuf` (in
  *Psychotria*'s). Their presence is itself a finding about plant mitogenomes, and they are
  exactly the MTPT trap that makes "a BLAST hit is a gene" wrong.
* `ccmFN2` (Arabidopsis' second copy) and the named ORFs `orfB` · `orfX`.

They are **kept in the reference set** with `-` in the vocabulary field: that way an annotation run
can *say* it found `rps12` instead of being silent about a gene it has no name for. Extending
`mito_gene_names.yaml` to cover the ten standard genes is a separate task — it changes naming for
every mito record and needs its own tests.

## Rebuilding

```
plastid panel pin --accessions NC_023209.1,NC_023103.1,NC_049064.1,NC_006581.1,NC_037304.1,\
NC_011033.1,NC_066984.1,NC_057654.1 --out .panel/mito-reference
```

then `build_proteins` / `build_rrna` / `render_fasta` from `plastid.core.mito_reference`. The check
that the committed FASTA still *is* what the panel builds is a test:
`tests/unit/test_mito_reference.py::test_the_committed_set_is_what_the_pinned_panel_builds`.
