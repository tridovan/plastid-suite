"""Reference Panel — propose comparison genomes in four taxonomy tiers (PROJECT.md §13.1).

Policy (decision §10.8): the engine only **proposes**, from a dated local cache of RefSeq
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

__all__ = [
    "IndexSnapshot",
    "PanelCandidate",
    "PanelEntry",
    "PanelProposal",
    "PinnedGenome",
    "load_index",
    "pin",
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
    accession: str
    path: str
    sha256: str
    md5: str
    bytes: int


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
        pinned.append(
            PinnedGenome(
                accession=entry.accession,
                path=str(path),
                sha256=hashlib.sha256(payload).hexdigest(),
                md5=hashlib.md5(payload).hexdigest(),
                bytes=len(payload),
            )
        )
    manifest = {
        "pinned_at": timestamp or datetime.now().isoformat(timespec="seconds"),
        "index_retrieved_at": retrieved_at,
        "genomes": [asdict(record) for record in pinned],
    }
    (directory / "PINNED.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return tuple(pinned)


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
