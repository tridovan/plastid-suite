# Release notes

What each version does, what it changes in a number you already have, and what it does not do. The
steps that cut a release are in [release.md](release.md); reproduced defects, each with the test that
holds it, are in [known issues](known_issues.md).

## 1.0.0 — not yet released

The first public version. Nothing was published before it, so there is nothing to upgrade from. The
migration table below is for anyone who ran an earlier checkout, because some numbers moved.

### What it does

- **A local interface.** `plastid ui` opens the toolkit in a browser on your own computer: records and
  the IR gate, analyses, assembly, references, results, figures and export. It starts the same
  command line you can type, shows that line beside each form, and sends nothing anywhere.
- **Reads to a checked sequence.** `plastid assemble` takes paired short reads (fastp, GetOrganelle)
  or PacBio HiFi reads (Oatk), maps the reads back with an aligner the assembler does not use, reads
  the graph for every structure that closes, and declares a sequence only when the reads cross every
  join and the graph is resolved. It runs here, or on a service of your choice through `prepare` and
  `import`; it never sends reads anywhere by itself.
- **What is there.** `plastid modules`, `plastid tools` and `plastid runs` list the analyses with their
  parameters, the tools this machine has, and the runs made, in text or JSON. A run is on disk from
  its first moment, and stopping one stops its tools and leaves a directory that says so.
- **Sequence to annotated record.** `plastid annotate` transfers from a directory of GenBank records
  you supply, after auditing them.
- **Records to tables.** Twenty-four modules ([the list, with what each costs](platforms.md)), among
  them structure and the IR gate, gene tables, codon usage, SSR, repeats, substitutions, diversity,
  Ka/Ks, population genetics, synteny, maximum-likelihood phylogeny and predicted RNA-editing sites.
- **Tables to what you hand in.** Figures, Word and Excel reports, a Supplementary package, a QC
  report and NCBI submission files.

Every run writes its own directory and records the name, version, command line and exit code of each
tool it called. An exported number carries its scope, and a value that was not measured says so.

### What it was checked against

| Path | Input | Result |
|---|---|---|
| Analyses | the golden pair *Iris domestica* NC_050833.1 vs *I. gatesii* NC_024936.1; 80 RefSeq plastomes in 10 families; 30 re-annotations; rotated records | the golden numbers in the [README](../README.md); module contracts in [module_contracts.md](module_contracts.md) |
| Short reads | six of the 50 read sets of Jin et al. 2020, each with a published answer, and *Arabidopsis thaliana* Col-0 | *Magnolia* and *Selaginella* identical to the record; *Arabidopsis* one base from the Sanger record; the published *Laurus* record has two places no read crosses ([getorganelle.md](tools/getorganelle.md)) |
| Assembly graph | 58 deposited graphs of that study | three cases with 12–15 structures are refused as unresolved; on *Picea abies* the published record is path 3 of 15 |
| HiFi reads | *Arabidopsis thaliana* Col-0 against the Sanger record AP000423.1 | closed at 154,480 bp, 2 bases from the record, 0 mismatches ([oatk.md](tools/oatk.md)) |
| RNA editing | *Vitis accD*, *Aegilops tauschii* KJ614412.1, *Nicotiana sylvestris* | *Vitis* reproduced exactly; *Aegilops* 32 of 34; *N. sylvestris* 31 of 33 with no false positive ([prep_cp.md](tools/prep_cp.md)) |

A passing test is not proof for another genome: the golden pair is one pinned input, which is why the
panels exist. The panels are third-party data and are not in the repository, so a clone can run
`make test` and `make golden` but not `make acceptance`.

CGAS 1.0.1 and CPStools 3.0 were run as comparisons. No module calls them and neither is a standard;
what each comparison shows, in both directions, is in [cgas.md](tools/cgas.md) and
[cpstools.md](tools/cpstools.md).

### If you ran an earlier checkout

| You ran | What changes now |
|---|---|
| A position-based module on a record with one SNP in a copy, an origin inside a copy, or no repeat annotation (K-01, K-03) | IR length, identity and region boundaries are measured on the full copies. A de-novo search used to return a fragment at 100 % |
| A sequence with `N` in the IR (K-04) | An `N` is no longer a match, so identity can fall and the gate can refuse |
| A bare FASTA with a one-base indel in an IR copy (K-07) | The detector used to report half the repeat and pass the gate. Re-run anything position-based |
| `diversity` or `substitutions` with a reference whose origin is not the LSC's first base (K-10) | The reference is normalized before it is aligned. Measured on the golden pair rotated by half: 386 substitutions became 1,094 |
| A cached run (K-02) | Older cache entries are not read. A moved or upgraded tool invalidates an entry |
| `codon_usage` or `amino_acid` on a record refused at the IR gate | They read only CDS sequences and now run, refusing by name a `/codon_start` other than 1 or a `/transl_table` other than 1 or 11 |
| Records spelling a name `tRNA-Ile`, `trnN (GUU)`, `rrn23s` or `tRNA-H (GUG)` (K-13) | They fold to the right locus. `tRNA-Ile` used to become `trnA-Ile`, an alanine tRNA |
| `gene_table` | `rpoC` is in the RNA polymerase group, and an `Unclassified` row is emitted and counts copies, so the rows sum to the total |
| `plastid panel pin` into a panel that already held records (K-16) | It used to replace the manifest. Check `PINNED.json` of any panel pinned twice |
| A script that reads exit codes (K-18) | An unknown module name now exits 2. A missing reference set and a missing fastp in `assemble prepare` now exit 3 |
| The R tree figure (K-06) | No tip is dropped and no outgroup is chosen unless the run names them; the default is a midpoint root |
| An outgroup or `drop-tips` for the R tree figure (K-28) | A name that is not a tip stops the figure and lists the tips; it used to be passed over. One outgroup tip roots the tree, and a named outgroup is the root without `--set root=outgroup` |
| The R figures of a run | They name genomes by species, from the records. The tables keep accessions, and the caption says which is which |
| `run --figures` on a run with `rna_editing` (K-24) | The RSCU, amino-acid and SSR figures are drawn. They used to be lost |
| A script that reads the exit code of a run with figures (K-25) | A figure that could be fed and was not drawn exits 3. It used to exit 0 |
| `phylogeny --figures` (K-26) | The run draws its own tree, without `--set tree=` |
| `--figures` with no plotting extra (K-29) | The analyses stay `ok`, and the error is in the run's `figures` section. It used to add a failed entry for the last analysis |
| A script that reads exit codes or `manifest.status` (K-23) | An input that cannot be read now exits 2, not 1. The status is `refused`, `tool_error` or `conflict` where it read `failed`, and a cancelled run exits 130 |
| A script that reads `manifest.command` (K-21) | It now holds what was typed, whole |
| `plastid run -m plastid_annotate` or `plastid batch -m plastid_annotate` (K-20) | Refused by name. Use `plastid annotate`, with `--set reannotate=yes` for a record that carries an annotation |
| `--set` with a misspelt word for a yes/no value (K-22) | It used to read as no. It is refused now |
| `--reference-set` naming one record instead of its folder (K-33) | Refused as a wrong command line, exit 2, naming the folder to pass. `annotate` and `reference audit` used to exit 3 and say the toolkit ships no records |
| `--exclude` with a file name or a comma list (K-34) | A file name leaves its record out; it used to stop the run. A record left out by a comma list is now left out of the manifest's `reference_sets` too |
| `annotate` with a reference file not named by its accession (K-35) | The manifest no longer counts the sample's own record as searched. A set whose every record is left out stops before the run starts, exit 2 |

The command line was reviewed before the first release, and names that meant two things were changed
without keeping the old ones, because nothing had been released under them:

| Command | Was | Now |
|---|---|---|
| `run`, `batch`, `annotate` | `--force` | `--no-cache`: it skips the stored results, and no check |
| `annotate` | `--reference DIR` | `--reference-set DIR` |
| `annotate` | `--audit-only`, `--no-confirm` | `plastid reference audit` |
| `reference verify`, `audit` | `--panel` | `--reference-set` |
| `reference loo` | `--panel`, `--out` | `reference leave-one-out`, `--reference-set`, `--to` |
| `assemble run` | `--reference` | `--compare-to` |
| `assemble run` | `--backend import` and its flags | `plastid assemble import` |
| `assemble run` | `--backend local` | `--skip-resource-check` |
| `ir-junction` | `--reference NAME` | `--order-from NAME` |
| `editing-reference build` | `--out FILE`, `--runs DIR` | `--to FILE`, `--out DIR` |
| `editing-reference verify` | `--reference` | `--reference-file` |
| `panel refresh` | `--out DIR`, `--retmax` | `--index FILE`, `--max-records` |
| `panel pin` | `--out` | `--to` |
| `panel propose` | `--mode A`, `--mode B` | `--mode all`, `--mode tiers` |
| `repeats` | `--set repeats=DIR` | `--set mummer=DIR` |

`--out` is now always the folder that holds runs. `plastid modules`, `plastid tools`, `plastid runs`
and `plastid ui` are new, and `annotate` takes `--label` like every other command that makes a run.

### What it does not do

The full list, with the reason for each, is in [things_we_do_not_do.md](things_we_do_not_do.md). The
ones a new user meets first:

- **Windows** runs through Docker only, and **macOS** has never been tested.
- **The Docker image** is built on your own computer and carries every tool but CPJSdraw, whose source
  states no licence: for IR junctions you name your own patched copy.
- **Nanopore reads** have no assembly route.
- **`rna_editing` is a prediction**, labelled PREP-Cp-compatible. It is not output of the original
  PREP-Cp, which ended in 2022, and it does not call editing from RNA-seq.
- **Oatk does not polish.** A one-base difference in a long homopolymer from the HiFi route is not a
  finding.
- **The toolkit ships no plastid records.** A reference set is yours to supply or to pin.
- **The tree is unrooted**, and a topology is not a verifiable number.
- **The mitochondrial path** is a frozen prototype.
