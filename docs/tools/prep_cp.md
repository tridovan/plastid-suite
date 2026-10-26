# PREP-Cp — selected RNA-editing predictor

**Decision, 26 October 2026:** use the original PREP-Cp for M4.5. PTS will prepare its inputs,
import its original outputs and map sites back to the sample. PTS will not implement a substitute
predictor or silently fall back to a conservation algorithm of its own. This is an integration
specification; no RNA-editing module has been implemented or exercised in this session.

## Selected scope and sources

The first integration reports **predicted plastid C-to-U sites in supported protein-coding genes**.
Input is an annotated GenBank in vocabulary A; RNA-seq is not required to run prediction. The
published PREP suite accepts spliced CDS, gene identity and starting codon position, supports batch
input and provides downloadable results/alignment files. See [Mower 2009, the original
method](https://doi.org/10.1093/nar/gkp337).

Use a declared cutoff of **0.8**, the value used in the selected public comparison study below.
Record it per query; never change it to force a target site count. Label all output `predicted`;
neither the score nor a site count is a measured editing frequency.

**Reference B** is PREP-Cp's own reference/alignment collection. PTS does not replace that
collection with its annotation DB or naming reference. PTS owns the mapping from sample locus/name
A to the predictor's supported gene identifiers, and from returned CDS positions to the original
sample coordinates. A user-defined alignment method would be a different integration and is not
selected here.

## Execution boundary

The original publication documents a web service. During this review the official
`prep.unl.edu` pages could not be accessed through the browsing tool. This does not prove the
service is permanently unavailable. No official standalone release, executable or supported CLI
has been verified; do not invent an installation command or assert a local version is available.

The selected integration therefore has two explicit stages around an external PREP-Cp run:

1. **Prepare:** write a batch submission and a query-to-sample manifest.
2. **Import:** read the results downloaded from the original service and the matching manifest.

The user can execute the external batch and retain its downloads. Automatic web submission is not
part of the initial adapter and must not be implemented against an assumed API. If an official
local distribution is subsequently obtained and verified, it may execute the same original engine
through `run_tool()`, with version/source provenance; that is not permission to reimplement it.

If the service cannot be used, retain the prepared input and report **awaiting external result**.
An existing genuine PREP-Cp output can be imported for the matching input; PTS must not synthesize
predictions, upload through another service or mark this run complete without a result. Parser
tests can use pinned genuine outputs, but end-to-end execution remains unverified until exercised.

## Input and output contracts

These are planned contracts, not commands that already exist.

| Artifact | Required contents |
|---|---|
| Input batch | Unique query ID, accepted gene identifier, spliced CDS, codon-start information and cutoff in the format documented by PREP-Cp |
| Query manifest | Sample/accession, raw GenBank hash, feature/locus identity, name in A, predictor gene ID, CDS hash, parts in transcript order, strand and coding-base-to-genome mapping |
| External-run record | Service URL, submission/retrieval date, parameters, exact input hash, available version/reference metadata and original result-file hashes |
| Raw outputs | Unmodified downloaded table, alignments and other supplied files; keep their original bytes |
| Site table | Sample/locus/name A, original tool gene label, CDS/genome positions, strand, reported edit/codon/amino-acid effect and raw score, with `predicted` scope |
| Coverage/status table | Every submitted, unsupported, ambiguous, rejected, absent-result or successfully evaluated query; zero predictions only when the tool actually evaluated that query |

If the service does not report its software or internal-reference revision, record **not reported**.
Do not substitute the article's year for a verified executable version or pretend the reference DB
has been pinned. Archived results permit repeatable import; they do not guarantee a future web run
will reproduce identical predictions.

Preparation must use transcript strand/order and declared frame. Import must validate the query,
CDS position, source base and supplied effect before mapping back. Minus-strand C-to-U corresponds
to a G-to-A change in the forward DNA frame; retain both conventions explicitly. Joined,
trans-spliced and origin-crossing features need controlled mapping tests. Ambiguous/unsupported
features are reported, not guessed. IR copies retain traceable feature locations; gene-based
prediction must not require a successful IR gate merely to extract an annotated CDS.

PTS never changes the input genome, `/translation` or submission annotation merely because a
site was predicted. Such a change would require a separate, explicit annotation decision.

## Selected public acceptance case

Use the *Aegilops tauschii* study by Wang et al., [DOI
10.3390/genes8010013](https://doi.org/10.3390/genes8010013). It specifies plastome **KJ614412.1**, a
PREP-Cp cutoff of 0.8, public RNA-seq experiments **SRX209402–SRX209425**, and site-level prediction,
RNA-read and experimental-validation results. Experiment accessions must be resolved to run/file
accessions when pinning; checksums and raw records have not been downloaded in this session.

Keep the paper's predictions, read-supported candidates and experimentally checked sites as
separate evidence sets. Compare sites in a common coordinate/gene frame. Unknown or unassessed
sites are not automatically negatives. Differences from the published totals need explanation;
they are not a reason to change thresholds or manufacture an identical count.

## Definition of Done

- Preparation and import consume the original engine's documented files; a captured genuine batch
  and output are pinned with provenance. No substitute algorithm is present.
- Imported raw scores/sites match the original output exactly; mapping tests cover both strands,
  joined/trans-spliced features, wraparound, duplicate gene names, stale input and missing queries.
- Unsupported/not-evaluated is distinguishable from zero predicted sites. Changed input or result
  files invalidate reuse; both the source run and the import are traceable in the manifest.
- The selected study is pinned and compared site by site with the correct evidence labels. Report
  agreement, false positives/negatives only where labels support them, and remaining unknowns.
- CLI/exports/reports use the same integration, carry `predicted` labels and retain raw downloads.
  Missing service/results and malformed imports have named statuses/errors.
- Exercise the original service or a verified original local distribution for end-to-end
  acceptance. Without that, only preparation/import can be accepted; M4.5 remains partly open.

RNA-seq calling with minimap2/samtools is **not selected as the M4.5 engine**. Its earlier simulation
is background evidence, not a replacement for PREP-Cp. RNA-seq from the selected study supplies an
independent check; supporting raw-read editing calls as a product feature would need its own scope.
