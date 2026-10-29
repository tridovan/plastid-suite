"""M11.10 — the reference flow from a fresh set, each step reading the previous one's artifact.

The DoD's chain, run as one test because the links are the claim:

1. **pin** two records from this tree with ``plastid panel pin --file``;
2. **audit** the pinned set (H-21 fixed: nothing blocks);
3. **choose** one record as the naming reference (vocabulary A);
4. **normalize** the other into it (``gene_normalize mode=normalize``) — the artifact must carry
   the sequence, locations, translations and the qualifiers outside the rename, and a new digest;
5. **annotate** the golden from the derived DB, re-read the artifact;
6. **verify** that artifact against the same naming reference — nothing may change.

The naming reference is deliberately the *original* record: it is already vocabulary A by
definition, and normalizing it against itself would prove nothing. No maintainer path is used:
the records come from ``tests/data``, the directory from this run's temporary tree.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

from conftest import require_tools
from plastid.cli import main
from plastid.core.genbank import read_genbank

GOLDEN = "NC_050833.1"
OTHER = "NC_024936.1"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _module_dir(run_root: Path, name: str) -> Path:
    found = sorted(run_root.glob(f"**/{name}"))
    assert found, f"{name} not found under {run_root}"
    return found[-1]


def _rows(module_dir: Path, table: str) -> list[dict[str, str]]:
    path = module_dir / "tables" / f"{table}.csv"
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    return list(csv.DictReader(line for line in lines if not line.startswith("#")))


@pytest.fixture(scope="module")
def flow(
    tools: dict[str, str | None],
    blast_tools: dict[str, str | None],
    data_dir: Path,
    tmp_path_factory: pytest.TempPathFactory,
) -> dict[str, Path]:
    require_tools(tools, "mafft")
    require_tools(blast_tools, "blastx", "blastn", "blastp", "makeblastdb", "aragorn")
    root = tmp_path_factory.mktemp("reference_flow")
    database = root / "reference"

    # 1 — pin the two records from this tree; the pin copies the bytes and records their digests
    assert (
        main(
            [
                "panel",
                "pin",
                "--file",
                str(data_dir / f"{GOLDEN}.gb"),
                "--file",
                str(data_dir / f"{OTHER}.gb"),
                "--out",
                str(database),
                "--note",
                "M11.10 flow fixture, copied from tests/data",
            ]
        )
        == 0
    )
    assert (database / "PINNED.json").is_file()

    # 2 — audit the pinned set: H-21 fixed, so the set has no blocked record
    audit_out = root / "audit"
    assert (
        main(
            [
                "reference",
                "audit",
                "--panel",
                str(database),
                "--no-confirm",
                "-o",
                str(audit_out),
            ]
        )
        == 0
    )

    # 3/4 — choose the naming reference (A) and normalize the other record into it
    normalize_out = root / "normalize"
    assert (
        main(
            [
                "run",
                "-m",
                "gene_normalize",
                "-i",
                str(database / f"{OTHER}.gb"),
                "--set",
                f"naming_reference={database / f'{GOLDEN}.gb'}",
                "-o",
                str(normalize_out),
            ]
        )
        == 0
    )
    normalize_dir = _module_dir(normalize_out, "01_gene_normalize")
    normalized = normalize_dir / f"{OTHER}.normalized.gb"
    assert normalized.is_file()

    # The derived DB: the naming reference as it came + the normalized record in A
    derived = root / "derived"
    derived.mkdir()
    (derived / f"{GOLDEN}.gb").write_bytes((database / f"{GOLDEN}.gb").read_bytes())
    (derived / f"{OTHER}.gb").write_bytes(normalized.read_bytes())

    # 5 — annotate the golden from the derived DB
    annotate_out = root / "annotate"
    assert (
        main(
            [
                "annotate",
                "-i",
                str(data_dir / f"{GOLDEN}.gb"),
                "--reference",
                str(derived),
                "--set",
                "reannotate=yes",
                "-o",
                str(annotate_out),
            ]
        )
        == 0
    )
    artifact = next(annotate_out.glob("**/plastid_annotation.gb"))

    # 6 — verify the artifact against the same naming reference: nothing may change
    verify_out = root / "verify"
    assert (
        main(
            [
                "run",
                "-m",
                "gene_normalize",
                "-i",
                str(artifact),
                "--set",
                f"naming_reference={database / f'{GOLDEN}.gb'}",
                "--set",
                "mode=verify",
                "-o",
                str(verify_out),
            ]
        )
        == 0
    )
    # 7 — the annotated artifact is the *input* of the next step, not a directory convention
    analysis_out = root / "analysis"
    assert main(["run", "-m", "genome_stats", "-i", str(artifact), "-o", str(analysis_out)]) == 0
    return {
        "database": database,
        "normalize_dir": normalize_dir,
        "normalized": normalized,
        "artifact": artifact,
        "verify_dir": _module_dir(verify_out, "01_gene_normalize"),
        "analysis_dir": _module_dir(analysis_out, "01_genome_stats"),
    }


def test_the_normalized_record_keeps_what_the_next_step_reads(flow, data_dir: Path) -> None:
    """Sequence, locations, translations and the qualifiers outside naming survive the rename."""
    source = read_genbank(data_dir / f"{OTHER}.gb")
    derived = read_genbank(flow["normalized"])

    assert derived.sequence == source.sequence
    # The name is what this step is allowed to change (vocabulary B → A); everything a later step
    # reads — kind, part order, coordinates, strand — must be the source's.
    assert [(g.kind, tuple(p.to_1based() for p in g.parts), g.strand) for g in derived.genes] == [
        (g.kind, tuple(p.to_1based() for p in g.parts), g.strand) for g in source.genes
    ]
    renamed = [
        (before.name, after.name)
        for before, after in zip(source.genes, derived.genes, strict=True)
        if before.name != after.name
    ]
    assert renamed, "the fixture was expected to exercise at least one rename"
    vocabulary_a = {gene.name for gene in read_genbank(data_dir / f"{GOLDEN}.gb").genes}
    assert all(after in vocabulary_a for _before, after in renamed), renamed
    assert [g.translation for g in derived.genes] == [g.translation for g in source.genes]
    for field in ("transl_table", "codon_start", "exception", "trans_splicing"):
        assert [getattr(g, field) for g in derived.genes] == [
            getattr(g, field) for g in source.genes
        ], field


def test_the_derived_record_has_a_new_digest_and_the_run_records_both(flow, data_dir: Path) -> None:
    original = _sha256(data_dir / f"{OTHER}.gb")
    normalized = _sha256(flow["normalized"])
    assert normalized != original, "a derived DB that is byte-identical was not converted"

    manifest = json.loads((flow["normalize_dir"].parent / "manifest.json").read_text())
    input_hashes = {entry["sha256"] for entry in manifest["inputs"]}
    file_hashes = {entry["sha256"] for entry in manifest["files"]}
    assert original in input_hashes
    assert normalized in file_hashes
    # the privileged record this step worked against is recorded by path and digest (D5-09)
    reference = manifest["modules"][0]["reference"]
    assert reference["accession"] == "NC_050833.1"
    assert reference["sha256"] == _sha256(flow["database"] / f"{GOLDEN}.gb")


def test_the_annotation_reads_back_and_carries_the_three_gene_classes(flow) -> None:
    annotated = read_genbank(flow["artifact"])
    counts = {
        kind: sum(1 for gene in annotated.genes if gene.kind == kind)
        for kind in ("CDS", "tRNA", "rRNA")
    }
    assert counts == {"CDS": 86, "tRNA": 38, "rRNA": 8}
    assert annotated.accession == GOLDEN


def test_verifying_the_annotated_artifact_changes_no_name(flow) -> None:
    summary = _rows(flow["verify_dir"], "normalize_summary")[0]
    assert summary["Mode"] == "verify"
    assert summary["Renamed"] == "0"
    assert summary["Conflicts"] == "0"


def test_the_annotated_artifact_is_the_input_of_the_next_step(flow) -> None:
    """M11.16's routing proof: the same bytes, named in both manifests.

    ``plastid run`` hands the modules the genomes it was given, so listing modules in one command
    is not a pipeline. The recipe chains the files instead: the annotated GenBank is the next
    run's input, and its digest is the same one the annotation run published as a file.
    """
    annotated_digest = _sha256(flow["artifact"])
    annotation_manifest = json.loads((flow["artifact"].parents[1] / "manifest.json").read_text())
    assert annotated_digest in {entry["sha256"] for entry in annotation_manifest["files"]}
    assert any(
        entry["path"].endswith("plastid_annotation.gb")
        for entry in annotation_manifest["artifacts"]
    )

    analysis_manifest = json.loads((flow["analysis_dir"].parent / "manifest.json").read_text())
    assert [entry["path"] for entry in analysis_manifest["inputs"]] == [str(flow["artifact"])]
    assert analysis_manifest["inputs"][0]["sha256"] == annotated_digest
