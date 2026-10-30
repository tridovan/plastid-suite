# The golden suite — what it proves, and its two tiers

Most numbers in plastome analysis have no true value to compare with, and this suite does not
pretend otherwise. Each check proves one of five things:

| Kind | Example | A pass shows | A pass does not show |
|---|---|---|---|
| A number measured once on the golden pair, then pinned | Ts/Tv 2.6917 (323 Ts / 120 Tv), π 0.007181, Ka/Ks 0.2167 | the number did not drift | that it is right — π and Ka/Ks have no published value for this pair |
| A second engine with the same definition | repeat-match and Vmatch report the same 9 repeats; the de-novo IR equals the annotated one | the result is not one implementation's bug | that the definition is the biologically right one |
| An answer known by construction, or an invariance | the synthetic IR record; five rotated origins; the reverse-complemented record; one sequence under two annotations | a property any correct tool must keep | — |
| A source outside this code | gene counts from the published record's own features; the curated mito record OX459128.1; NCBI's VADR tables | agreement with that source | that the source is right — a disagreement is investigated on both sides |
| The pinned tools really run | MAFFT strategy names, the MACSE jar, table2asn | the argv and versions the code claims | — |

The golden pair (*Iris domestica* NC_050833.1 vs *Iris gatesii* NC_024936.1) is one input. It is
not evidence about any other genome; the panel tier is.

## Two tiers

| Command | Selects | Needs | Cost | Runs |
|---|---|---|---|---|
| `make golden` | `tests/golden`, not marked `panel` | the pinned tools | 85 cases, 5 min 32 s | on every push in CI, and locally |
| `make acceptance` | every test marked `panel`, unit and golden | `.panel/` and the pinned tools | 234 cases, ~40 min | locally, when closing a milestone |

The 5 min 32 s was measured on the development machine (WSL2) on 29 October 2026, with every tool
present. The ~40 min is an estimate: it is the 43-minute full golden of the same day minus the
`make golden` share. The first `make acceptance` run prints its own durations.

A test is marked `panel` when it reads `.panel/`: 80 RefSeq plastomes in 10 families, 30 PGA
re-annotations, rotated FASTA and the pinned reference sets. The directory is gitignored, so a
fresh clone and CI cannot run that tier. `make acceptance` refuses to start without it and sets
`PTS_REQUIRE_TOOLS=1` and `PTS_REQUIRE_PANEL=1`, so a missing tool, extra or panel subset fails
instead of skipping. Both targets print the 25 slowest tests.

A new test that reads `.panel/` gets `@pytest.mark.panel`, or the mark in its module's
`pytestmark`. Without it, the test skips in CI and belongs to no tier that says so.

A changed golden number is never edited to make the suite pass: explain it with independent
evidence and get the dev's review first.
