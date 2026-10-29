"""The tool inventory, kept in step with the code and the three environment files (M11.19, G-04).

**The defect this guards.** `environment.yml` once mixed the tools an analysis runs with
assemblers, dating programs and engines that had only been surveyed. An installed package reads
like a capability: five PGA runs that produced zero genes would have looked like success to any
pipeline checking that an output file exists, and the same mistake at the environment level is
"the toolkit can assemble a plastome" inferred from `getorganelle=` in a dependency list.

**The claim that replaces it.** `docs/tools/README.md` gives every tool one role —
`active-runtime` · `cross-check` · `future` · `survey-only` — and a role decides which environment
file carries it. These checks hold that claim from four directions:

* an adapter in `src/plastid/external/` with no row fails;
* a `cross-check` row whose adapter some module imports fails — cross-check means *nothing calls
  it*, not "we also happen to have it";
* an `active-runtime` row with no call site fails, because that is a promise with nothing behind
  it;
* a conda package in any environment file that no row claims, or in the wrong file for its role,
  fails — both directions, so neither the table nor the files can drift alone.

Each helper is also exercised against fabricated input that must fail it: a check nobody has
watched fail is not a check.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "plastid"
README = ROOT / "docs" / "tools" / "README.md"

ENVIRONMENTS = {
    "active-runtime": ROOT / "environment.yml",
    "cross-check": ROOT / "environment-check.yml",
    "future": ROOT / "environment-future.yml",
    "survey-only": ROOT / "environment-future.yml",
}
ROLES = ("active-runtime", "cross-check", "future", "survey-only")

# A package a program ships inside another conda package (PAML's `mcmctree`) is not a separate
# line, so the package cell may only name what an environment file can actually install.
_CONDA_NAME = re.compile(r"^[a-z0-9][a-z0-9-]*[a-z0-9]$")
_ROW = re.compile(r"^\|\s*`([^`]+)`\s*\|([^|]*)\|([^|]*)\|([^|]*)\|([^|]*)\|")
_ADAPTER = re.compile(r"external/([a-z_]+)\.py")
_NONE = ("—", "-", "none", "")


def _inventory(text: str) -> dict[str, dict[str, str]]:
    """The tool rows of the inventory table, keyed by tool name."""
    try:
        section = text.split("## The inventory", 1)[1].split("\n## ", 1)[0]
    except IndexError:
        return {}
    rows: dict[str, dict[str, str]] = {}
    for line in section.splitlines():
        match = _ROW.match(line)
        if match is None:
            continue
        name, role, package, version, call_site = (part.strip() for part in match.groups())
        rows[name] = {
            "role": role,
            "package": package,
            "version": version,
            "call_site": call_site,
        }
    return rows


def _packages(cell: str) -> list[str]:
    """Conda package names a cell names; an em dash (with or without a note) names none.

    The table backticks every tool name, so a package cell reads `` `r-base`, `r-ape` `` — the
    backticks are markup and are stripped here rather than made part of the claim.
    """
    if not cell or cell in _NONE or cell.startswith("—"):
        return []
    names = (part.strip().strip("`") for part in cell.split(","))
    return [name for name in names if name and name not in _NONE]


def _conda_packages(path: Path) -> dict[str, str]:
    """``name -> spec`` for the conda dependencies of one environment file.

    Only the ``dependencies:`` block counts — ``channels:`` is a list of the same shape, and
    reading it as dependencies would make every channel look like an undeclared tool.
    """
    packages: dict[str, str] = {}
    inside = False
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.split("#", 1)[0].rstrip()
        if not stripped.strip():
            continue
        if not stripped.startswith((" ", "\t")) and stripped.endswith(":"):
            inside = stripped.strip() == "dependencies:"
            continue
        if not inside or not stripped.strip().startswith("- "):
            continue
        spec = stripped.strip()[2:].strip()
        name = re.split(r"[=<>!]", spec, maxsplit=1)[0].strip()
        if not _CONDA_NAME.match(name):
            continue
        packages[name] = spec
    packages.pop("python", None)  # the interpreter every file needs, not a declared tool
    return packages


def _problems(
    rows: dict[str, dict[str, str]], adapter_names: set[str], environments: dict[str, str]
) -> list[str]:
    """Every way the table and the code/environment files disagree, as readable lines."""
    problems: list[str] = []

    for name, row in rows.items():
        if row["role"] not in ROLES:
            problems.append(f"{name}: role {row['role']!r} is not one of {ROLES}")
        if row["role"] == "active-runtime":
            call_site = row["call_site"]
            if call_site in _NONE or call_site.startswith("—"):
                problems.append(
                    f"{name}: claims active-runtime but names no call site — either a module "
                    f"calls it and the row says where, or it is not runtime"
                )

    named = {adapter for row in rows.values() for adapter in _ADAPTER.findall(row["call_site"])}
    for adapter in sorted(adapter_names - named):
        problems.append(
            f"external/{adapter}.py is not in the inventory — a tool the code can call must "
            f"have a role (docs/tools/README.md)"
        )

    # package -> the role(s) whose rows claim it. Two roles may share a file (future and
    # survey-only both live in environment-future.yml), so membership is by role, not by file.
    claims: dict[str, set[str]] = {}
    for row in rows.values():
        for package in _packages(row["package"]):
            claims.setdefault(package, set()).add(row["role"])
    installed = {
        package for path in set(environments.values()) for package in _conda_packages(path)
    }

    for package in sorted(installed - set(claims)):
        problems.append(
            f"a conda environment installs {package!r}, which no inventory row claims — a "
            f"package in an environment is not a feature"
        )
    for package, roles in sorted(claims.items()):
        if package not in installed:
            if any(
                _conda_packages(environments[role]).keys() for role in roles if role in environments
            ):
                problems.append(
                    f"{package}: no environment file installs it, but a row claims it — the "
                    f"install line and the table disagree"
                )
            continue
        for role in sorted(roles):
            path = environments.get(role)
            if path is None or package in _conda_packages(path):
                continue
            problems.append(
                f"{package}: role {role} belongs in {path.name}, which does not install it — "
                f"the package sits in a file for a different role"
            )
    return problems


def _importers(adapter: str, src: Path = SRC) -> list[str]:
    """Product files that import the adapter — ``modules``, the CLI and the figure layer."""
    needle = f"plastid.external.{adapter}"
    hits: list[str] = []
    for path in sorted(src.rglob("*.py")):
        if path.parent.name == "external":
            continue
        if needle in path.read_text(encoding="utf-8"):
            hits.append(path.relative_to(src.parent).as_posix())
    return hits


def _cross_check_calls(rows: dict[str, dict[str, str]], src: Path = SRC) -> list[str]:
    """Adapters a `cross-check` row names that a product file nonetheless imports."""
    problems: list[str] = []
    for name, row in rows.items():
        if row["role"] != "cross-check":
            continue
        for adapter in _ADAPTER.findall(row["call_site"]):
            callers = _importers(adapter, src)
            if callers:
                problems.append(
                    f"{name} is cross-check, but external/{adapter}.py is imported by "
                    f"{', '.join(callers)} — cross-check means no module calls it"
                )
    return problems


# ------------------------------------------------------------------ the real repository


def test_every_adapter_and_environment_package_has_one_role() -> None:
    rows = _inventory(README.read_text(encoding="utf-8"))
    assert rows, "the inventory table moved and this check would pass over nothing"
    adapters = {path.stem for path in (SRC / "external").glob("*.py") if path.stem != "__init__"}
    assert adapters, "no adapters found — the path moved, and this check would pass vacuously"
    problems = _problems(rows, adapters, ENVIRONMENTS)
    assert not problems, "tool inventory disagrees with the code/environment:\n" + "\n".join(
        problems
    )


def test_no_cross_check_tool_is_imported_by_a_product_module() -> None:
    """The distinction G-04 is about, measured rather than asserted: Vmatch is the one engine the
    `repeats` module must never depend on, because its licence is academic-only."""
    rows = _inventory(README.read_text(encoding="utf-8"))
    problems = _cross_check_calls(rows)
    assert not problems, "\n".join(problems)


def test_every_active_runtime_adapter_is_actually_imported() -> None:
    """A runtime row names a call site; the call site must exist.

    The opposite direction of the cross-check guard: there, a `cross-check` adapter must have no
    importer; here, an `active-runtime` adapter must have one. Without this, a row could keep
    claiming a command path after the code stopped calling the tool — the table would describe a
    capability the package no longer has.
    """
    rows = _inventory(README.read_text(encoding="utf-8"))
    missing: list[str] = []
    for name, row in rows.items():
        if row["role"] != "active-runtime":
            continue
        for adapter in _ADAPTER.findall(row["call_site"]):
            if not _importers(adapter):
                missing.append(f"{name} (external/{adapter}.py)")
    assert not missing, "runtime tools whose adapter no product file imports: " + ", ".join(missing)


def test_the_runtime_check_reports_an_adapter_nobody_imports(tmp_path: Path) -> None:
    """The negative control: a runtime row for an adapter removed from the code must fail."""
    fake_src = tmp_path / "src"
    (fake_src / "modules").mkdir(parents=True)
    rows = {
        "ghost": {
            "role": "active-runtime",
            "package": "—",
            "version": "—",
            "call_site": "external/vmatch.py",
        }
    }
    missing = [
        adapter
        for row in rows.values()
        if row["role"] == "active-runtime"
        for adapter in _ADAPTER.findall(row["call_site"])
        if not _importers(adapter, fake_src)
    ]
    assert missing == ["vmatch"]


def test_every_runtime_conda_package_is_version_pinned() -> None:
    """CLAUDE.md: a tool a module uses is declared *with a pinned version*. An unpinned MAFFT
    would hand two users two Ts/Tv numbers from one genome.

    The R plotting packages are the one deliberate exception and the inventory says why in the
    row: they draw a picture and never produce a number that reaches a table, so a version bump
    cannot move a measured value. That exception is named here, not assumed — any *other*
    unpinned runtime package fails.
    """
    unpinned_by_design = {
        "r-ggplot2",
        "r-dplyr",
        "r-tidyr",
        "r-reshape2",
        "r-rcolorbrewer",
        "r-zoo",
        "r-scales",
        "r-pheatmap",
        "r-ape",
    }
    rows = _inventory(README.read_text(encoding="utf-8"))
    unpinned: list[str] = []
    runtime = _conda_packages(ROOT / "environment.yml")
    for name, row in rows.items():
        if row["role"] != "active-runtime":
            continue
        for package in _packages(row["package"]):
            if package in unpinned_by_design:
                continue
            spec = runtime.get(package)
            if spec is not None and "=" not in spec:
                unpinned.append(f"{name} ({package}: {spec})")
    assert not unpinned, "runtime tools installed without a version pin: " + ", ".join(unpinned)


def test_the_minimum_install_is_still_two_packages() -> None:
    """The extras gate themselves inside functions; the minimum install must not grow a tool."""
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    block = pyproject.split("dependencies = [", 1)[1].split("]", 1)[0]
    names = {line.strip().strip('",').split(">=")[0].lower() for line in block.splitlines()}
    names.discard("")
    assert names == {"biopython", "pyyaml"}, (
        f"the minimum install is {sorted(names)} — CLAUDE.md says biopython + PyYAML alone run "
        f"every pure-Python analysis"
    )


# ----------------------------------------------------------------------- negative controls


def test_the_inventory_check_reports_an_undocumented_adapter() -> None:
    rows = {
        "mafft": {
            "role": "active-runtime",
            "package": "mafft",
            "version": "7.526",
            "call_site": "external/mafft.py",
        }
    }
    problems = _problems(rows, {"mafft", "newtool"}, {})
    assert any("external/newtool.py is not in the inventory" in line for line in problems)


def test_the_inventory_check_reports_a_runtime_row_with_no_call_site() -> None:
    rows = {
        "ghost": {"role": "active-runtime", "package": "—", "version": "—", "call_site": "none"}
    }
    problems = _problems(rows, set(), {})
    assert any("claims active-runtime but names no call site" in line for line in problems)


def test_the_inventory_check_reports_an_unknown_role() -> None:
    rows = {"x": {"role": "maybe", "package": "—", "version": "—", "call_site": "—"}}
    problems = _problems(rows, set(), {})
    assert any("role 'maybe' is not one of" in line for line in problems)


def test_the_environment_check_reports_an_unclaimed_package(tmp_path: Path) -> None:
    environment = tmp_path / "environment.yml"
    environment.write_text(
        "dependencies:\n  - python>=3.10\n  - mafft=7.526\n  - sibelia=3.0.7\n",
        encoding="utf-8",
    )
    rows = {
        "mafft": {
            "role": "active-runtime",
            "package": "mafft",
            "version": "7.526",
            "call_site": "external/mafft.py",
        }
    }
    problems = _problems(rows, {"mafft"}, {"active-runtime": environment})
    assert any("installs 'sibelia'" in line for line in problems)


def test_the_environment_check_reports_a_package_in_the_wrong_file(tmp_path: Path) -> None:
    """A deferred assembler installed into the *runtime* file is exactly G-04's original shape."""
    runtime = tmp_path / "environment.yml"
    runtime.write_text("dependencies:\n  - mafft=7.526\n  - sibelia=3.0.7\n", encoding="utf-8")
    future = tmp_path / "environment-future.yml"
    future.write_text("dependencies:\n  - python>=3.10\n", encoding="utf-8")
    rows = {
        "mafft": {
            "role": "active-runtime",
            "package": "mafft",
            "version": "7.526",
            "call_site": "external/mafft.py",
        },
        "sibelia": {"role": "future", "package": "sibelia", "version": "3.0.7", "call_site": "—"},
    }
    problems = _problems(rows, {"mafft"}, {"active-runtime": runtime, "future": future})
    assert any("sibelia" in line and "different role" in line for line in problems), problems


def test_the_cross_check_check_reports_a_real_caller(tmp_path: Path) -> None:
    """Fabricate a product file that imports Vmatch — the check must notice."""
    fake_src = tmp_path / "src"
    (fake_src / "modules").mkdir(parents=True)
    (fake_src / "modules" / "repeats.py").write_text(
        "from plastid.external.vmatch import run_vmatch\n", encoding="utf-8"
    )
    rows = {
        "vmatch": {
            "role": "cross-check",
            "package": "vmatch",
            "version": "2.3.1",
            "call_site": "external/vmatch.py",
        }
    }
    problems = _cross_check_calls(rows, src=fake_src)
    assert any("cross-check means no module calls it" in line for line in problems)


@pytest.mark.parametrize("cell", ["— (obtained separately)", "none", ""])
def test_an_em_dash_package_cell_names_no_conda_package(cell: str) -> None:
    assert _packages(cell) == []
