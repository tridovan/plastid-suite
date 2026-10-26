"""The plastid reference set an annotation is transferred from (P3B-01, D-b).

A plastome is annotated by transferring a curated annotation from relatives — the method the
mito path already uses (``core/mito_reference.py``) and the one D5-08 settles for the plastid
path. This module owns the *set*: which records are in it, that it still is what was pinned,
and — for validation — which records may annotate which target.

**The set is chosen by the taxonomy rule** (D-b): same family first, the rule the Reference
Panel and the mito reference set already follow. The rule is applied when the set is *built*
(a person pins a family's records, with the pin recording what was taken), not per run: a run
takes 1..N reference records and the engine searches them.

**Leave-one-out is the validation design, written down here because it is a rule and not a
script.** The corpus (P4a) is 34 Iridaceae plastomes whose published
annotation is the answer key. Annotating a record with a reference set that contains that very
record would be the annotation grading its own homework — the corpus's PGA set carries the
golden's own sequence and reads 85 CDS where the contract says 86. So for every record
in the corpus, the reference set that annotates it is **the same set minus itself**, and
:func:`leave_one_out` is that rule as data: target → the records that may annotate it.

**What this module does not do.** It does not select references per run, does not read the
records' genes, and does not rename anything — the annotator (P3B-07) and the naming step
(D5-09, P3B-04) do those. Here a reference set is its pin, its rule, and its integrity.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from plastid.core.panel import load_pinned
from plastid.core.resources import data_file

__all__ = [
    "MANIFEST_FILENAME",
    "RECORDED_MANIFEST",
    "ReferenceSet",
    "build_manifest",
    "leave_one_out",
    "load_recorded_manifest",
    "read_manifest",
    "resolve_directory",
    "verify_pin",
    "write_manifest",
]

# The manifest lives beside PINNED.json, in the reference set's own directory: a pin says what
# the bytes are, the manifest says what the set is for and which target each record may annotate.
MANIFEST_FILENAME = "REFERENCE_SET.json"

# The *committed* copy of that manifest. The records themselves live outside the repository
# (`.panel/` is gitignored, and a 34-record plastome set is data, not source), so this file is what
# a run reads to learn where they are and which of them may annotate a given target.
RECORDED_MANIFEST = "plastid_reference.json"


@dataclass(frozen=True, slots=True)
class ReferenceSet:
    """A pinned set of reference records, the rule that chose them, and the LOO design.

    ``family`` is the taxonomy level the rule was applied at, and ``rule`` is the sentence a
    later reader needs to judge whether the set still answers the question it was built for.
    ``leave_one_out`` maps each target accession to the accessions allowed to annotate it.
    """

    directory: str
    family: str
    rule: str
    built_at: str
    records: tuple[str, ...]
    leave_one_out: dict[str, tuple[str, ...]]

    def references_for(self, target: str) -> tuple[str, ...]:
        """The records that may annotate ``target`` — its leave-one-out set, or all of them.

        A target outside the set is annotated by the whole set: nothing of it is inside, so
        nothing has to be left out.
        """
        return self.leave_one_out.get(target, self.records)

    @property
    def families(self) -> frozenset[str]:
        return frozenset({self.family})


def leave_one_out(accessions: Iterable[str]) -> dict[str, tuple[str, ...]]:
    """target → the same records minus the target, sorted, for every record in the set.

    A repeated accession is refused rather than de-duplicated: two files claiming one accession
    would make "minus itself" remove both, and the run would silently search a smaller set than
    the manifest promised.
    """
    ordered = sorted(accessions)
    seen: set[str] = set()
    repeated: list[str] = []
    for accession in ordered:
        if accession in seen:
            repeated.append(accession)
        seen.add(accession)
    if repeated:
        raise ValueError(
            f"the reference set carries {sorted(set(repeated))} more than once — "
            f"leave-one-out cannot say which copy to exclude"
        )
    return {
        accession: tuple(other for other in ordered if other != accession) for accession in ordered
    }


def verify_pin(directory: str | Path) -> tuple[str, ...]:
    """Re-hash every pinned record; one ``file — reason`` line per problem, and nothing else.

    Three ways a set stops being the set that was pinned, each reported by name: a file whose
    bytes no longer hash to what ``PINNED.json`` records, a file the pin names but that is gone,
    and a ``*.gb`` sitting in the directory that the pin never listed — a record nobody checked
    is a record no run may search (the same failure the corpus fixtures' SHA256SUMS guard
    catches, tests/golden/test_corpus.py).
    """
    root = Path(directory)
    records = load_pinned(root)
    problems: list[str] = []
    listed: set[str] = set()
    for record in records:
        listed.add(Path(record.path).name)
        path = _pinned_file(root, record.path)
        if path is None:
            problems.append(f"{Path(record.path).name} — in PINNED.json but not in the directory")
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != record.sha256:
            problems.append(
                f"{path.name} — sha256 {digest[:12]}… but the pin records {record.sha256[:12]}…"
            )
    for path in sorted(root.glob("*.gb")):
        if path.name not in listed:
            problems.append(f"{path.name} — in the directory but not in PINNED.json")
    return tuple(problems)


def _pinned_file(root: Path, recorded: str) -> Path | None:
    """Where a pinned record actually is, from the path the pin recorded.

    ``pin_files`` writes the path it was given, so a pin made with a relative ``--out`` records
    a path relative to *that* shell's working directory. The record's home is the panel either
    way — the recorded path is tried first so a pin that points somewhere else is honoured, and
    the panel itself is the fallback that keeps ``plastid reference verify`` working from any
    directory. Neither resolving means the file is gone, which the caller reports.
    """
    path = Path(recorded)
    for candidate in (path, root / path.name):
        if candidate.is_file():
            return candidate
    return None


def build_manifest(
    directory: str | Path,
    *,
    family: str,
    rule: str,
    built_at: str | None = None,
) -> ReferenceSet:
    """Read a pinned set and write down what it is: the family, the rule, and the LOO design."""
    root = Path(directory)
    records = tuple(record.accession for record in load_pinned(root))
    return ReferenceSet(
        directory=str(root),
        family=family,
        rule=rule,
        built_at=built_at or datetime.now().isoformat(timespec="seconds"),
        records=records,
        leave_one_out=leave_one_out(records),
    )


def write_manifest(reference_set: ReferenceSet, path: str | Path) -> Path:
    """Write the manifest where a run finds it — text, sorted, one accession per line."""
    payload = {
        "directory": reference_set.directory,
        "family": reference_set.family,
        "rule": reference_set.rule,
        "built_at": reference_set.built_at,
        "records": list(reference_set.records),
        "leave_one_out": {
            target: list(references)
            for target, references in sorted(reference_set.leave_one_out.items())
        },
    }
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    return target


def load_recorded_manifest() -> ReferenceSet | None:
    """The committed copy of the set's design, or ``None`` where this tree has none.

    A tree without the recorded manifest is not a broken tree — the minimum install ships no
    reference set, and the annotate step then says so by name rather than guessing a directory.
    """
    path = data_file(RECORDED_MANIFEST)
    return read_manifest(path) if path.is_file() else None


def resolve_directory(given: str | Path | None) -> Path:
    """Where the reference records are: what the run was told, or the recorded set's own directory.

    One source of truth for the default, because two callers need it and they must not disagree:
    the annotator (which searches the records) and ``plastid annotate --audit-only`` (which audits
    them without annotating anything) — a flag that audited a *different* directory than the run
    would be worse than no flag.
    """
    if given is not None:
        return Path(given)
    recorded = load_recorded_manifest()
    return Path(recorded.directory if recorded is not None else "")


def read_manifest(path: str | Path) -> ReferenceSet:
    """Read a committed manifest back — the offline half of the set, usable without the pin."""
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return ReferenceSet(
        directory=str(payload["directory"]),
        family=str(payload["family"]),
        rule=str(payload["rule"]),
        built_at=str(payload["built_at"]),
        records=tuple(payload["records"]),
        leave_one_out={
            str(target): tuple(references)
            for target, references in payload["leave_one_out"].items()
        },
    )
