# The plastid reference set — what it is, why these records, and the leave-one-out rule

**Why a set exists at all.** An assembled plastome arrives as a bare sequence (the assembler's
output), so the annotate step transfers an annotation from curated relatives — the method the
mito path already uses (`src/plastid/data/mito_reference.md`) and the one D5-08 settles for the
plastid path. This is the set that transfer searches.

## The rule that chose the records (D-b)

**Same family first** — the rule the Reference Panel and the mito reference set already follow.
For the plastid path the family is **Iridaceae**, applied at pin time (a set is *built* for a
family; a run takes the records it is given). Two reasons this is the family, both of them about
what this project already has to check against:

- the **contract genome is Iridaceae** — *Iris domestica* `NC_050833.1` and its comparison genome
  *Iris gatesii* `NC_024936.1` (PROJECT.md §6);
- the **answer key is Iridaceae** — the Điềm corpus is 34 Iridaceae plastomes with a **published
  annotation** on every one, which is what makes a transfer result checkable at all (§8.31).

A sample from another family is annotated by a set pinned for *its* family, by the same rule —
the set is per-family data, not a per-run choice.

## Sources — pinned 2026-10-08, NCBI RefSeq plastomes, complete and annotated

The 34 records are the pinned Iridaceae corpus (`.panel/iridaceae-34/`, sha256 for sha256,
PROJECT.md §8.31), re-pinned as the reference set: `.panel/plastid-reference/`, **34/34 passing
the IR gate**, 5,204,491 bp in all. Every sha256 below was measured on the file in that
directory, not copied from a table.

| Accession | Organism | Length (bp) | CDS | sha256 |
|---|---|---|---|---|
| `NC_024936.1` | *Iris gatesii* | 153,441 | 86 | `fd4417273a23…` |
| `NC_041459.1` | *Crocus cartwrightianus* | 150,819 | 86 | `17f36ee80972…` |
| `NC_041460.1` | *Crocus sativus* | 150,820 | 86 | `f113a0aeb0f2…` |
| `NC_042827.1` | *Iris missouriensis* | 153,084 | 86 | `bd8f9b34b8e5…` |
| `NC_050833.1` | *Iris domestica* | 153,724 | 86 | `a946f23e6e5f…` |
| `NC_062351.1` | *Gladiolus x gandavensis* | 151,654 | 84 | `c27ee3d64cd6…` |
| `NC_072561.1` | *Aristea ecklonii* | 153,548 | 86 | `33e43e6ca199…` |
| `NC_072562.1` | *Chasmanthe floribunda* | 151,293 | 86 | `abc14e204ddf…` |
| `NC_072563.1` | *Phalocallis coelestis* | 152,789 | 86 | `ae7b281f2096…` |
| `NC_072564.1` | *Dierama igneum* | 152,328 | 86 | `ea23428a423b…` |
| `NC_072565.1` | *Diplarrena moraea* | 155,895 | 86 | `b16107d5d186…` |
| `NC_072566.1` | *Freesia laxa* | 152,858 | 85 | `94f458fbf5ef…` |
| `NC_072572.1` | *Gladiolus communis* | 151,608 | 86 | `982604e2c637…` |
| `NC_072573.1` | *Trimezia candida* | 151,102 | 86 | `55e8f56a87f4…` |
| `NC_072574.1` | *Isophysis tasmanica* | 158,654 | 86 | `f13110ed201c…` |
| `NC_072575.1` | *Libertia pulchella* | 154,129 | 87 | `2e2c35b88d96…` |
| `NC_072576.1` | *Moraea polystachya* | 152,361 | 86 | `0a30c4339a4b…` |
| `NC_072577.1` | *Moraea spathulata* | 152,127 | 86 | `474c85040816…` |
| `NC_072578.1` | *Nivenia stokoei* | 152,780 | 86 | `0112e9bd5cbd…` |
| `NC_072579.1` | *Olsynium douglasii* | 152,197 | 86 | `02c004a624df…` |
| `NC_072580.1` | *Patersonia fragilis* | 155,062 | 86 | `f47070088d38…` |
| `NC_072581.1` | *Sisyrinchium idahoense* | 151,006 | 86 | `15a03059ee77…` |
| `NC_072582.1` | *Tigridia pavonia* | 164,622 | 86 | `4f9cac33cdb2…` |
| `NC_072583.1` | *Tritonia disticha* | 151,724 | 86 | `9ec0f6b916dd…` |
| `NC_072584.1` | *Watsonia pillansii* | 151,981 | 86 | `34941360d4d6…` |
| `NC_072585.1` | *Witsenia maura* | 151,260 | 86 | `2835c6ee43ae…` |
| `NC_072661.1` | *Hesperantha coccinea* | 151,955 | 86 | `8924698cf536…` |
| `NC_072662.1` | *Trimezia steyermarkii* | 152,777 | 86 | `2c836e148ffc…` |
| `NC_081929.1` | *Iris afghanica* | 153,234 | 86 | `fc4495a5ec0d…` |
| `NC_081930.1` | *Iris bismarckiana* | 153,187 | 86 | `e99a8d559039…` |
| `NC_081931.1` | *Iris hermona* | 153,185 | 86 | `24d0e2d18d18…` |
| `NC_081932.1` | *Iris lortetii* | 153,194 | 86 | `bf131479d905…` |
| `NC_088562.1` | *Iris orchioides* | 151,622 | 87 | `4d0ab83fae0c…` |
| `PX516262.1` | *Eleutherine bulbosa* | 152,471 | 87 | `a7f3baa0c006…` |

`PX516262.1` is a GenBank submission, not a RefSeq record — it is the tool-paper species and the
corpus carries it, so it is in the set. The **CDS column ranges 84–87**: 132 = 86 + 38 + 8 is the
golden genome's count, not a constant, and a reference set whose members disagree about a gene's
existence is the point of check 3 below, not a defect to hide.

**One caveat, stated because a manuscript may not carry it silently.** These records are the
corpus — the project's *validation* dataset (D-b: "the corpus is validation, never production").
They are the reference pool here because the contract genome is Iridaceae and no other Iridaceae
source is pinned; a sample of another family gets a set pinned for that family by the same rule.
**No number measured against this set becomes a manuscript table for a real sample.**

## The leave-one-out rule (the validation design)

For a transfer result to mean anything, the reference set that annotates a record must not
contain that record: the corpus carries the golden's own sequence annotated by PGA, and PGA's
reading is 85 CDS where the contract says 86 (§8.31) — the answer key inside the reference set
would let an annotator agree with itself.

So, for every record in the set, the reference set that annotates it is **the same 34 minus
itself** — 33 records. Machine-readable in `plastid_reference.json` (committed, so the design is
readable without the 28 MB pin): `records` is the set; `leave_one_out[target]` is what may
annotate that target. A target *outside* the set is annotated by the whole set — nothing of it is
inside, so nothing has to be left out.

## Rebuilding

```
# 1. pin the records (source: the corpus already in the tree)
plastid panel pin $(for f in .panel/iridaceae-34/*.gb; do printf -- "--file %s " "$f"; done) \
  --out .panel/plastid-reference --note "…"
# 2. re-verify every sha256 — `file — reason` per problem, exit 1 when the set changed
plastid reference verify --panel .panel/plastid-reference
# 3. write the leave-one-out manifest (and the committed offline copy)
plastid reference loo --panel .panel/plastid-reference --family Iridaceae \
  --out src/plastid/data/plastid_reference.json
```

The checks that the committed manifest still describes the pin, and that every sha256 of the pin
matches, are tests: `tests/unit/test_plastid_reference.py`.
