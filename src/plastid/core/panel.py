"""Reference Panel — propose comparison genomes in four taxonomy tiers.

Policy: the engine only **proposes**, from a dated local cache of RefSeq
records. The user picks the final list; pinning downloads the chosen genomes and records
their checksums, so every later run works offline and stays reproducible.

Tiers: 0 same species · 1 same genus · 2 same family (different genus) · 3 same order
(different family). Mode A takes everything at one level (family by default); mode B takes
up to ``per_tier`` representatives per tier.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Iterable
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

from plastid.core.genbank import read_genbank
from plastid.core.ir import IRDetectionError, IRGateError, detect_ir, enforce_ir_gate
from plastid.core.model import Genome, count_genes

__all__ = [
    "IndexSnapshot",
    "PanelCandidate",
    "PanelEntry",
    "PanelProposal",
    "PinnedGenome",
    "entries_from_summaries",
    "load_index",
    "pin",
    "pin_files",
    "propose",
    "save_index",
    "tier_of",
]

TIER_LABELS = {
    0: "same species",
    1: "same genus",
    2: "same family, different genus",
    3: "same order, different family",
}

DEFAULT_MIN_LENGTH = 100_000  # a complete plastome; shorter records are fragments


@dataclass(frozen=True, slots=True)
class PanelEntry:
    accession: str
    organism: str
    species: str = ""
    genus: str = ""
    family: str = ""
    order: str = ""
    length: int = 0
    annotated: bool = False
    refseq: bool = True


@dataclass(frozen=True, slots=True)
class IndexSnapshot:
    retrieved_at: str
    entries: tuple[PanelEntry, ...]


@dataclass(frozen=True, slots=True)
class PanelCandidate:
    entry: PanelEntry
    tier: int

    @property
    def tier_label(self) -> str:
        return TIER_LABELS.get(self.tier, f"tier {self.tier}")


@dataclass(frozen=True, slots=True)
class PanelProposal:
    target: PanelEntry
    mode: str  # "A" | "B"
    level: str  # the level used by mode A
    candidates: tuple[PanelCandidate, ...]
    notes: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class PinnedGenome:
    """One pinned download, with everything measured by reading the file back.

    The measurement fields are *data*: a diverged or absent IR is recorded, never enforced
    here — the run-time gate stays the only enforcement point (code seam 4).
    """

    accession: str
    path: str
    sha256: str
    md5: str
    bytes: int
    # where the bytes came from: "ncbi" (fetched by accession) or "file" (produced elsewhere —
    # a web service, a cluster, a colleague's machine). A heavy tool runs where it fits
    # (P1-12); either way its answer arrives as a file, and a file that is not
    # pinned is a number nobody can reproduce.
    source: str = "ncbi"
    note: str | None = None  # which service/when, for a pinned file
    length: int | None = None  # read back from the file — None when it failed to parse
    cds: int | None = None
    trna: int | None = None
    rrna: int | None = None
    topology: str | None = None
    ir_length: int | None = None
    ir_source: str | None = None  # "annotation" | "denovo_diag"; None when no IR was found
    ir_identity: float | None = None
    ir_error: str | None = None  # IRDetectionError (absent) or IRGateError (diverged)
    parse_error: str | None = None


def entry_from_summary(summary: dict, lineage: str) -> PanelEntry:
    """Build a PanelEntry from an NCBI esummary record plus its taxonomy lineage.

    An NCBI lineage for a species ends "...; <order>; <family>; <genus>" — the last component
    is the genus, the one before it the family, and the one before that the order.
    """
    parts = [part.strip() for part in lineage.split(";") if part.strip()]
    genus = parts[-1] if parts else ""
    family = parts[-2] if len(parts) >= 2 else ""
    order = parts[-3] if len(parts) >= 3 else ""
    organism = str(summary.get("organism", "")).strip()
    accession = str(summary.get("accessionversion") or summary.get("caption") or "").strip()
    if not accession:
        raise ValueError(f"esummary record has no accession: {summary!r}")
    refseq = accession.startswith("NC_")
    return PanelEntry(
        accession=accession,
        organism=organism,
        species=organism,
        genus=genus or (organism.split(" ")[0] if organism else ""),
        family=family,
        order=order,
        length=int(summary.get("slen") or 0),
        # RefSeq is annotated by NCBI's build. A GenBank submission's annotation is only known
        # after download (measured at pin time), so it is never assumed here — the old
        # unconditional True was a claim nothing had checked.
        annotated=refseq,
        refseq=refseq,
    )


def entries_from_summaries(
    records: Iterable[tuple[dict, str]],
) -> tuple[tuple[PanelEntry, ...], tuple[str, ...]]:
    """``(esummary record, lineage)`` pairs → entries, plus one problem line per unreadable record.

    Unreadable records are returned, never dropped silently — the caller prints them and fails
    loudly (CLAUDE.md: fail loud). Shared by ``panel refresh`` and ``panel pin``.
    """
    entries: list[PanelEntry] = []
    problems: list[str] = []
    for summary, lineage in records:
        try:
            entries.append(entry_from_summary(summary, lineage))
        except ValueError as error:
            problems.append(str(error))
    return tuple(entries), tuple(problems)


def tier_of(target: PanelEntry, entry: PanelEntry) -> int | None:
    """The tier an entry belongs to relative to the target; None when unrelated or itself."""
    if entry.accession == target.accession:
        return None
    if entry.species and entry.species == target.species:
        return 0
    if entry.genus and entry.genus == target.genus:
        return 1
    if entry.family and entry.family == target.family:
        return 2
    if entry.order and entry.order == target.order:
        return 3
    return None


def _passes_quality(entry: PanelEntry, min_length: int) -> bool:
    """The three manual criteria of mode A: RefSeq · complete genome · annotated."""
    return entry.refseq and entry.length >= min_length and entry.annotated


def propose(
    snapshot: IndexSnapshot,
    target_accession: str,
    *,
    mode: str = "B",
    level: str = "family",
    per_tier: int = 3,
    min_length: int = DEFAULT_MIN_LENGTH,
) -> PanelProposal:
    """Propose comparison genomes for one target — the user still picks the final list."""
    target = _find(snapshot.entries, target_accession)
    notes = [f"index retrieved {snapshot.retrieved_at}; engine proposes, the user confirms"]
    scored: list[PanelCandidate] = []
    dropped = 0
    for entry in snapshot.entries:
        tier = tier_of(target, entry)
        if tier is None:
            continue
        if not _passes_quality(entry, min_length):
            dropped += 1
            continue
        scored.append(PanelCandidate(entry=entry, tier=tier))
    if dropped:
        notes.append(f"{dropped} candidate(s) dropped: not RefSeq, too short, or unannotated")

    if mode == "A":
        ceiling = {"species": 0, "genus": 1, "family": 2, "order": 3}.get(level)
        if ceiling is None:
            raise ValueError(f"mode A needs level in species/genus/family/order, got {level!r}")
        # "Everything at family level" includes same-genus and same-species records.
        candidates = [candidate for candidate in scored if candidate.tier <= ceiling]
        notes.append(
            f"mode A: everything at {level} level or closer ({len(candidates)} candidate(s))"
        )
    elif mode == "B":
        candidates = []
        for tier in sorted({candidate.tier for candidate in scored}):
            in_tier = sorted(
                (candidate for candidate in scored if candidate.tier == tier),
                key=lambda candidate: candidate.entry.accession,
            )
            candidates.extend(in_tier[:per_tier])
        notes.append(f"mode B: up to {per_tier} representatives per tier")
    else:
        raise ValueError(f"unknown mode {mode!r} — use 'A' or 'B'")

    candidates.sort(key=lambda candidate: (candidate.tier, candidate.entry.accession))
    return PanelProposal(
        target=target, mode=mode, level=level, candidates=tuple(candidates), notes=tuple(notes)
    )


def _measure_ir(genome: Genome) -> tuple[int | None, str | None, float | None, str | None]:
    """``(length, source, identity, error)`` for a freshly read genome — the gate as data.

    A mitogenome legitimately has no IR and a bad assembly legitimately fails the gate; both
    are *recorded* here and only enforced later, at run time (code seam 4).
    """
    try:
        ir = detect_ir(genome)
    except IRDetectionError as error:
        return None, None, None, f"{type(error).__name__}: {error}"
    try:
        enforce_ir_gate(ir)
    except IRGateError as error:
        return ir.length, ir.source, ir.identity, f"{type(error).__name__}: {error}"
    return ir.length, ir.source, ir.identity, None


def pin(
    entries: Iterable[PanelEntry],
    fetch: Callable[[str], str],
    out_dir: str | Path,
    *,
    timestamp: str | None = None,
    retrieved_at: str | None = None,
) -> tuple[PinnedGenome, ...]:
    """Download the chosen genomes once, record checksums, and index them in ``PINNED.json``.

    ``fetch(accession) -> GenBank text`` is injected so tests never touch the network.
    """
    directory = Path(out_dir)
    directory.mkdir(parents=True, exist_ok=True)
    pinned: list[PinnedGenome] = []
    for entry in entries:
        text = fetch(entry.accession)
        payload = text.encode("utf-8")
        path = directory / f"{entry.accession}.gb"
        path.write_bytes(payload)
        pinned.append(_measured(entry.accession, path, payload))
    manifest = {
        "pinned_at": timestamp or datetime.now().isoformat(timespec="seconds"),
        "index_retrieved_at": retrieved_at,
        "genomes": [asdict(record) for record in pinned],
    }
    (directory / "PINNED.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return tuple(pinned)


def _measured(
    accession: str,
    path: Path,
    payload: bytes,
    *,
    source: str = "ncbi",
    note: str | None = None,
) -> PinnedGenome:
    """One pinned record, with everything a later run needs measured by reading the file back.

    A file that does not re-read is still pinned — with ``parse_error`` set rather than an
    exception: refusing to record it would leave a number in a report whose source nobody kept.
    """
    length: int | None = None
    cds: int | None = None
    trna: int | None = None
    rrna: int | None = None
    topology: str | None = None
    ir_length: int | None = None
    ir_source: str | None = None
    ir_identity: float | None = None
    ir_error: str | None = None
    parse_error: str | None = None
    try:  # verify the pin by reading it back with our own reader
        genome = read_genbank(path)
        counts = count_genes(genome)
        length = len(genome)
        cds, trna, rrna = counts.cds, counts.trna, counts.rrna
        topology = genome.topology
        ir_length, ir_source, ir_identity, ir_error = _measure_ir(genome)
    except Exception as error:  # noqa: BLE001 — recorded, not swallowed
        parse_error = f"{type(error).__name__}: {error}"
    return PinnedGenome(
        accession=accession,
        path=str(path),
        sha256=hashlib.sha256(payload).hexdigest(),
        md5=hashlib.md5(payload).hexdigest(),
        bytes=len(payload),
        source=source,
        note=note,
        length=length,
        cds=cds,
        trna=trna,
        rrna=rrna,
        topology=topology,
        ir_length=ir_length,
        ir_source=ir_source,
        ir_identity=ir_identity,
        ir_error=ir_error,
        parse_error=parse_error,
    )


def pin_files(
    paths: Iterable[str | Path],
    out_dir: str | Path,
    *,
    note: str | None = None,
    timestamp: str | None = None,
) -> tuple[PinnedGenome, ...]:
    """Pin records **produced somewhere else** — a web service, a cluster, another machine.

    This is the other half of the plan's rule for a heavy tool (P1-12): it runs
    where it fits, on this machine when it can and elsewhere when it cannot, and its answer
    arrives as a *file*. Pinning that file is what makes the answer reproducible — the same
    bytes, sha256 and md5 as an NCBI download, the same re-read with our own reader, and the
    same IR measurements.

    Two differences from :func:`pin`: there is no accession to fetch, so the key is the record's
    own id — or the file's stem when it carries none, which is exactly the case for a sample that
    has no accession yet — and an existing ``PINNED.json`` is **added to**, not replaced, because
    a file is normally pinned *into* a panel that already exists.
    """
    directory = Path(out_dir)
    directory.mkdir(parents=True, exist_ok=True)
    pinned: list[PinnedGenome] = []
    for raw in paths:
        source_path = Path(raw)
        payload = source_path.read_bytes()
        accession = _record_id(source_path)
        target = directory / f"{accession}.gb"
        if source_path.resolve() != target.resolve():
            target.write_bytes(payload)
        pinned.append(_measured(accession, target, payload, source="file", note=note))

    existing = _load_pinned(directory)
    keys = {record.accession for record in pinned}
    merged = [record for record in existing if record.accession not in keys]
    manifest = {
        "pinned_at": timestamp or datetime.now().isoformat(timespec="seconds"),
        "index_retrieved_at": None,
        "genomes": [asdict(record) for record in (*merged, *pinned)],
    }
    (directory / "PINNED.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return tuple(pinned)


def _record_id(path: Path) -> str:
    """What to call a pinned file: the record's own id, else the file's name.

    The record's id wins because that is what every other tool in the pipeline keys on; the file
    name is the fallback for an annotation of a sample that has no accession yet.
    """
    try:
        genome = read_genbank(path)
    except Exception:  # noqa: BLE001 — an unreadable file still gets pinned, and reported
        return path.stem
    return genome.accession or path.stem


def _load_pinned(directory: Path) -> tuple[PinnedGenome, ...]:
    """The records already in a directory's PINNED.json — empty when there is none."""
    manifest = directory / "PINNED.json"
    if not manifest.is_file():
        return ()
    try:
        payload = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return ()
    known = {field for field in PinnedGenome.__dataclass_fields__}
    return tuple(
        PinnedGenome(**{key: value for key, value in record.items() if key in known})
        for record in payload.get("genomes", [])
        if isinstance(record, dict) and "accession" in record
    )


def save_index(path: str | Path, snapshot: IndexSnapshot) -> None:
    payload = {
        "retrieved_at": snapshot.retrieved_at,
        "entries": [asdict(entry) for entry in snapshot.entries],
    }
    Path(path).write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def load_index(path: str | Path) -> IndexSnapshot:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    entries = tuple(PanelEntry(**entry) for entry in payload.get("entries", []))
    return IndexSnapshot(retrieved_at=str(payload.get("retrieved_at", "")), entries=entries)


def load_pinned(out_dir: str | Path) -> tuple[PinnedGenome, ...]:
    """Read ``PINNED.json`` — the offline path used by every later run."""
    payload = json.loads((Path(out_dir) / "PINNED.json").read_text(encoding="utf-8"))
    return tuple(PinnedGenome(**record) for record in payload.get("genomes", []))


def _find(entries: tuple[PanelEntry, ...], accession: str) -> PanelEntry:
    for entry in entries:
        if entry.accession == accession:
            return entry
    raise KeyError(f"{accession!r} is not in the panel index — refresh the cache first")
