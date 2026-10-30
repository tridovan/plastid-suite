# PREP-Cp-compatible RNA-editing prediction

**Decision, 29 October 2026.** PTS runs the published PREP method itself and reports
**PREP-Cp-compatible** predicted C-to-U sites in supported plastid protein-coding genes. A site is
never presented as output of the original PREP-Cp, because the original can no longer be run (see
below). The original reference alignments are not being requested from the author.

**Implemented 30 October 2026** as the module `rna_editing` (`plastid run -m rna_editing`) and the
command `plastid editing-reference build | verify`. It reproduces both archived genuine PREP-Cp
outputs exactly, reproduces 32 of the 34 sites the original predicted for *Aegilops tauschii*, and
finds 31 of the 33 known sites of the unseen *Nicotiana sylvestris* with no prediction outside them.
Every difference is recorded below with the reference column that explains it.

## Why the original PREP-Cp cannot be used

Checked on 29 October 2026:

| Finding | Evidence |
|---|---|
| The service ended | The PREP home page, Internet Archive snapshot of 19 April 2022: "The PREP suite web service will be ending June 30, 2022." It points users to Deepred-Mt, a **mitochondrial** predictor |
| The host is gone | `prep.unl.edu` has no DNS record, while `unl.edu` resolves |
| No local release existed | Both papers give only the web address ([Mower 2005](https://doi.org/10.1186/1471-2105-6-96), [Mower 2009](https://doi.org/10.1093/nar/gkp337)), and the archived site has no download page |
| The reference was never published | The 2009 paper has no supplementary data. The 2005 paper published PREP-Mt's 42 mitochondrial alignments; the 35 chloroplast alignments existed only on the server, and the Internet Archive holds no copy |
| Independent confirmation | [He et al. 2026](https://doi.org/10.3389/fpls.2026.1916637) state that the server and its manually adjusted alignments are no longer available, and report a "PREP-Cp-compatible" reconstruction whose scripts are available only on request |

Other plastid predictors were checked and not chosen:

- PREPACT 3 is web-only, with no local release or API, and uses its own reference-based method.
- CURE-Chloroplast was web-only, and its server no longer resolves.
- PlantC2U publishes code without a licence and uses a different method (a convolutional network).

Deepred-Mt and PREP-Mt are mitochondrial predictors and must not be applied to plastid genes.

## The method

The algorithm is from Mower 2005; the chloroplast reference and the cutoff are from Mower 2009. The
rule is `src/plastid/core/rna_editing.py`; the module only wires the aligner and writes the tables.

1. **Input:** the spliced CDS of one supported gene, read from `/codon_start`, and a cutoff C
   between 0 and 1.
2. **Translate** with the standard genetic code. Its amino-acid assignments are those of plastid
   table 11; only the start-codon sets differ. The terminal stop is removed.
3. **Align** the translation to the gene's reference protein alignment, without re-aligning the
   reference sequences among themselves. PREP used a ClustalW profile alignment; PTS adds the query
   with MAFFT 7.526 through `run_tool()`:
   `mafft --amino --localpair --maxiterate 1000 --quiet --add query.fasta reference.fasta`.
   The reference rows must come back unchanged once the columns the query inserted are removed;
   this is checked on every call, and a locus that fails it is `not evaluated` (alignment failure).
4. **Score** every codon that contains a C. List the unedited state and every C-to-U edited state
   (for `CCG`: `CCG`, `UCG`, `CUG`, `UUG`). A state's score is the fraction of reference residues
   in that alignment column, gaps excluded, that equal the state's amino acid. Scores are exact
   fractions; the matching and compared counts are exported beside every score.
5. **Choose a state.** The highest score wins. On a tie, the state with fewer edits wins, so a
   silent edit is never predicted. Between equally edited states, the one edited at codon
   position 2 wins.
6. **Report** a site only when the chosen state is edited and its score is **≥ C**. Mower 2005
   defines the threshold as "greater than or equal to". The 2009 paper's "> C" is a looser
   wording, and the archived genuine output below reports a 0.80 site at C = 0.8. C is compared as
   a fraction parsed from its decimal text: `0.8` read as a binary float is slightly more than 4/5,
   which would drop exactly that site.

**Supported genes.** These are the 35 genes PREP-Cp offered: accD, atpA, atpB, atpF, atpI, ccsA,
clpP, matK, ndhA, ndhB, ndhD, ndhF, ndhG, petB, petD, petG, petL, psaB, psaI, psbB, psbE, psbF,
psbL, rpl2, rpl20, rpl23, rpoA, rpoB, rpoC1, rpoC2, rps2, rps8, rps14, rps16 and ycf3. A CDS
counts as one of them by its name in vocabulary A — its canonical name or any synonym of its locus
(`pafI` reaches `ycf3`). Any other gene is reported as **not supported**, never as zero sites.

**Cutoff 0.8.** Mower 2009 measured 0.8 as the optimum for seed-plant chloroplasts: MCC 76 % in
leave-one-out tests on the training taxa, and 82 % on unseen data. The cutoff is recorded per run
(`--set cutoff=…`) and is never changed to reach a target site count. Every output is
`predicted`: neither a score nor a site count is a measured editing frequency.

**The open points, decided.**

- **Gaps do not count.** Within one gene, the archived C = 0.2 output reports scores in fifths
  (0.60, 0.80) next to scores in quarters (0.75, 0.50). The denominator therefore changes by
  column: the original counted only the reference sequences with a residue in that column. Mower
  2005's formula leaves gaps unstated, and He et al. 2026 made the same choice.
- **A tie between single edits at codon positions 1 and 3 cannot change a call.** Under the
  standard code a C-to-U change at codon position 3 never changes the amino acid (NNC and NNU
  always encode one amino acid), so any state that carries a position-3 edit ties with the same
  state without it — which has fewer edits and wins. The code still orders the remaining states
  totally, by their edited positions; `tests/unit/test_rna_editing.py` checks the argument on all
  64 codons.
- **The MAFFT argv** is the line in step 3, chosen as the accuracy-oriented L-INS-i option for at
  most eight sequences, and used for building B′ too. It is in every table's scope label.

A codon containing an IUPAC ambiguity code is not evaluated, and is listed as such; the rest of the
locus is scored. A C codon whose column holds no reference residue (the query's own insertion) is
listed as unscored.

## Reference B′ — rebuilt by PTS, never called the original reference

**Taxa.** The seven PREP-Cp training taxa, one record each, pinned by SHA-256 in
`src/plastid/data/rna_editing/taxa.yaml`. The files are kept out of git:
`plastid panel pin --accessions <the seven> --out .panel/rna_editing` re-creates them byte for byte,
and the build refuses any other bytes. Mower 2009
names the last two only by genus; the species are the ones whose plastid editing had been surveyed
by then. These taxa belong to the method; they are not a toolkit default for any other analysis.

| Taxon | Record | Edit sites in B′ | In PREP-Cp's own reference (Mower 2009, Table 2: TP + FN) | Edit-list source |
|---|---|---|---|---|
| *Arabidopsis thaliana* | NC_000932.1 | 31 | 29 | Chateigner-Boutin & Small 2007, Table 1 — the complete 2007 list |
| *Nicotiana tabacum* | NC_001879.2 | 33 | 38 | Sasaki et al. 2003, Table 2 — includes the 31 sites of Hirose et al. 1999 |
| *Atropa belladonna* | NC_004561.1 | 30 | 31 | Schmitz-Linneweber et al. 2002 — 28 of tobacco's 31 sites edited in *Atropa*, plus 3 of its own |
| *Zea mays* | NC_001666.2 | 24 | 26 | Wang et al. 2017, Table 4; Chen et al. 2017 Supp. Table 4 (rpoC2) |
| *Oryza sativa* (Japonica) | NC_001320.1 | 23 | 25 | Wang et al. 2017, Table 4 |
| *Phalaenopsis aphrodite* subsp. *formosana* | NC_007499.1 | 41 | 42 | Zeng et al. 2007, as tabulated with genome positions by Chen et al. 2017, Supp. Table 2 |
| *Pinus thunbergii* | NC_001631.1 | 27 | 27 | Wakasugi et al. 1996, Table 1 (26); accD P88L from Sasaki Y. et al. 2001 through UniProtKB/Swiss-Prot |

**How the lists compare with 2009's.** Five of the seven come from the surveys PREP-Cp was trained
on, read in their articles; the maize and rice lists come from a later comparison table (Wang et
al. 2017), because Maier et al. 1995 and Corneille et al. 2000 were not read. The counts differ from
Mower 2009's in both directions and no site is inferred from another species. A missing edit leaves
the unedited residue in that taxon's row, which lowers the edited state's score in that column: it
can cost a prediction, never add one — as it does at the rpoC2 site below, where rice and *Pinus*
have no rpoC2 edit sourced. Changing a list is a data change: rows in `edit_sites.tsv` with their
source, a rebuild, a new digest.

**Every site is sourced and checked.** `edit_sites.tsv` gives each site its source, the place in
the source and the source's own coordinate. Every site was placed on its pinned record and checked
there: the base is a C on the transcript strand, and the codon and both amino acids are the ones
the source states. A source whose own fields disagree is kept only where the record settles it,
and the row says how (a genome coordinate one base off its own codon fields; a CDS model that
starts elsewhere, so the codon number differs while the codon is the same; an amino acid the
printed codon cannot give). Four such cases settle on the records: tobacco rpoC2 is numbered on
the 1986 model, which starts at the in-frame ATG four codons upstream of the RefSeq one (in Z00044.1
codon 1248 is the TCA); black pine ORF320 is ccsA (ten of its eleven listed codons are the record's
ccsA codons); black pine petD's stop-creating CAA is codon 161, as the paper's 160-residue product
requires, not the printed 160; and *Atropa*'s three "absent" sites are T in its record, as the
paper says. Sites outside the 35 genes, non-coding sites and silent edits are not part of B′.
Three supported CDS that their records name only by `/product` (*Phalaenopsis* psaI and petL,
*Pinus* petL) are named in `taxa.yaml` by location and product, with the reason; *Pinus* petL's
model is corrected there too — the record starts it one codon after the ACG that editing makes its
start, and its own `/transl_except` puts Met on a CCT.

**Build.** For each supported gene, extract the CDS of each taxon that carries it (the first copy
when there are two), apply that taxon's sourced edits, translate with the standard code, drop the
terminal stop, and align the mature proteins once with MAFFT 7.526 L-INS-i. Each taxon contributes
one protein per gene: accD has five rows (the two grasses lack it) and rpoC1 seven, which is the
fifths and sevenths of the archived scores. A stop created by a sourced edit ends the mature
protein (black pine petD: 160 residues, as in the other taxa); a taxon whose protein would keep an
internal stop no source created is excluded for that gene, with the reason recorded (none is,
today).

**One digest.** It is the SHA-256 of the canonical JSON of everything B′ is made of — record pins,
edit sites, sources, recipe, aligner version and argv, and the aligned rows. It is in every table's
scope label, in the run manifest (`reference_sets`, role `rna_editing_reference`) and in the step's
parameters, so the cache key changes when B′ does. The packaged file is
`src/plastid/data/rna_editing/reference_b_prime.json`; a file edited after its build is refused.

```sh
plastid panel pin --accessions NC_000932.1,NC_001879.2,NC_004561.1,NC_001666.2,NC_001320.1,NC_007499.1,NC_001631.1 --out .panel/rna_editing
plastid editing-reference verify --records .panel/rna_editing
plastid editing-reference build --out <file>
```

`build` refuses a record whose bytes are not its pin and a MAFFT that is not 7.526, and writes its
own run (manifest and tool logs). Rebuilding today reproduces the packaged file byte for byte
(digest `a2162d886cdfba70b8e68842736cdc05385ccac77fab456569874a65edf558b9`).

**Expected differences from the original.** A different aligner, no manual adjustment, and edit
lists that differ from 2009's in both directions. Differences are measured and explained, never
tuned away. If the original alignments ever become public, they enter as a separate reference with
their own label. Results from the two references are never reported under one label.

## Contracts

| Artifact | Contents |
|---|---|
| `rna_editing_query` | Sample, sequence SHA-256 (the raw file's hash is in the run manifest's `inputs`), name in the record and in A, supported gene, CDS hash, parts in transcript order, strand, `/codon_start`, `/transl_table`, partial flag, and the coding-base-to-genome mapping as segments |
| `rna_editing_sites` | Sample, locus and name A; CDS, amino-acid and alignment positions; genome position (1-based), strand, the base and the change on the forward strand; codon and amino acid before and after; score with its matching and compared counts; reference digest; cutoff; status `predicted · PREP-Cp-compatible` |
| `rna_editing_coverage` | Every CDS and every supported gene with one status: evaluated (with its site count, which may be 0), not supported, absent from the record, or not evaluated with the reason (pseudogene, translation table, frame, internal stop, alignment failure); ambiguous and unscored codons listed |
| `alignments/<locus>.fasta` | The query-plus-reference alignment of each evaluated locus, so any score can be recomputed by hand |

Mapping rules:

- The CDS is read the way `coding_sequence` reads it: every part on its own strand, in annotation
  order. Codon and amino-acid numbers start at `/codon_start`; CDS positions count from the
  feature's first base.
- A minus-strand C-to-U appears as G-to-A on the forward strand; both conventions are exported.
- Joined and origin-crossing features map through their parts (tested on hand-built records).
- A gene with copies in both IRs is evaluated and mapped once per copy (`ndhB#1`, `ndhB#2`); two
  features that share a name are two loci too. Prediction is gene-based and does not need the IR
  gate.
- PTS never changes the input genome, `/translation` or the submission annotation because a site
  was predicted.

## Acceptance evidence

All checks run in `tests/golden/test_rna_editing_acceptance.py` (tier B, MAFFT 7.526); the
rebuild and the *N. sylvestris* test read `.panel/` and are marked `panel`.

1. **Archived genuine outputs (regression) — reproduced exactly.** The input is the sample on the
   archived PREP-Cp form: the *Vitis vinifera* accD CDS of NC_007957.1 (1,527 nt, identical base
   for base), first nucleotide at codon position 1. The Internet Archive kept two original PREP-Cp
   results for it, captured between January and March 2018:

   | Nt Pos | AA Pos | Align Col | Effect | Score at C = 0.2 | Reported at C = 0.8 | PTS score (column in B′) |
   |---|---|---|---|---|---|---|
   | 125 | 42 | 42 | ACG (T) => ATG (M) | 1.00 | yes | 1.00 = 4/4 (42) |
   | 134 | 45 | 45 | CCT (P) => CTT (L) | 0.75 | — | 0.75 = 3/4 (45) |
   | 185 | 62 | 62 | ACG (T) => ATG (M) | 0.50 | — | 0.50 = 1/2 (62) |
   | 244 | 82 | 85 | CAT (H) => TAT (Y) | 0.75 | — | 0.75 = 3/4 (85) |
   | 484 | 162 | 168 | CAC (H) => TAC (Y) | 0.75 | — | 0.75 = 3/4 (169) |
   | 845 | 282 | 311 | TCG (S) => TTG (L) | 0.80 | yes | 0.80 = 4/5 (313) |
   | 992 | 331 | 372 | ACA (T) => ATA (I) | 0.60 | — | 0.60 = 3/5 (382) |
   | 1436 | 479 | 520 | CCA (P) => CTA (L) | 0.60 | — | 0.60 = 3/5 (530) |

   PTS reports these eight sites and no other at C = 0.2, and exactly the two at C = 0.8, with the
   same positions, effects and scores. The 0.80 at 845 needs two sourced edits, *Arabidopsis*
   accD S265L and *Pinus* accD P88L; 0.60 at 1436 needs *Phalaenopsis* accD P471L. The alignment
   column depends on the reference and is not compared. A third archived output, for rpoC1, came
   from an input that was not kept, so it cannot serve as a regression.

2. ***Aegilops tauschii* KJ614412.1** ([Wang et al.](https://doi.org/10.3390/genes8010013)),
   C = 0.8, three evidence sets kept apart (`tests/data/rna_editing/wang2017_aegilops_sets.tsv`,
   each row placed on the record and checked). PTS predicts 34 sites in 15 genes (each IR copy
   is mapped; a site is counted once per gene).

   - **Table 2, the original PREP-Cp prediction (34 sites in 15 genes): 32 reproduced.** Not
     reproduced, each with the column that decides it: rpoC2 4031 (L 5/7 — the S are rice and
     *Pinus*, which have no rpoC2 edit sourced: rice's comparison table has no rpoC2 row and the
     1996 black-pine survey did not examine rpoC2), and rpoC2 2158 (S 1/2), a grass-only insertion
     in B′'s alignment with five gap rows, which the unpublished hand-adjusted original cannot be
     compared on. Predicted by PTS but not in Table 2: petB 20 and rpoC2 2123, both
     single-residue columns (maize's petB is annotated 19 codons longer at its N-terminus; the
     rpoC2 site is in the grass-only insertion).
   - **RNA-seq (Table 3, 60 calls):** 26 are sites the method can predict (one of the 35 genes, an
     amino-acid change, not a stop); PTS predicts 22 of them. Not predicted: rpoC2 4031 (above),
     and ndhG 347, petL 56 and rpoA 527, which the original server did not predict either.
   - **RT-PCR (24 sites):** 23 confirmed edited, 20 of them predictable; PTS predicts 17, as the
     original did (both miss ndhG 347, petL 56 and rpoA 527). `rpoB`-617, which RT-PCR showed
     unedited, is predicted at 0.86 by PTS and was predicted by the original: a conserved-residue
     false positive of the method, not of this reference.

   A site nobody assessed is not a negative, and "not called by RNA-seq" is reported as that, not
   as a false prediction.

3. **Mower 2009 Table 3 — *Nicotiana sylvestris*, unseen.** The 33 sites of Sasaki et al. 2003's
   *N. sylvestris* column are exactly the 33 Mower scored (his TP 32 + FN 1); placed on NC_007500.1
   (`tests/data/rna_editing/sasaki2003_nsylvestris_sites.tsv`). At C = 0.8 PTS predicts 31 of them
   and nothing else (Mower: TP 32, FN 1, FP 2). The two it misses: rps2 134 (I 4/7 — the T are
   *Atropa*, never assessed at this 2003 site by the 2002 survey, and the two grasses, whose source
   has no rps2 row) and rpoC2 3731 (the *Aegilops* rpoC2 column above). *Pisum* and *Saccharum*:
   not measured — their edit lists were not sourced.

RNA-seq is not an input to prediction; it serves only as the independent check above. Calling
editing from reads would be a separate feature with its own scope. The earlier minimap2/samtools
simulation is background evidence, not this engine.
