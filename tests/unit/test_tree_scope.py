"""A phylogeny module would need three things this toolkit does not have (M11.21, proposal).

**Why this file exists before any code.** A tree cannot be drawn from nothing: `r_tree` draws a
`<file>.treefile` the user supplies, and no module produces one. Before promising a module, the
inputs it would consume have to exist, and today two of the three do not.

**The three requirements, measured 28/10/2026:**

1. **A reusable alignment.** `diversity`, `popgen` and `kaks` each build an `Alignment` in memory
   through `external/mafft.align()` and export **tables**, never the alignment itself: a grep for
   `"fasta"`/`".aln"` under `src/plastid/modules/` and `src/plastid/result/` returns nothing. A
   tree step therefore cannot read what those modules already computed; it would re-align from
   scratch, which is a second numeric path and exactly what M11.16's artifact contract forbids.
2. **A concatenated, per-gene supermatrix.** The paper's Methods use *PhyloSuite v1.2.2* "screen
   the PCGs" — i.e. a partitioned PCG supermatrix — and this toolkit has no concatenation code
   (`concatenat`/`supermatrix`/`partition` have no product hits). `popgen` aligns per gene and
   reports per-gene statistics; it does not join them into one matrix.
3. **An inference tool, declared and pinned.** `.dev` records that CGAS's module 14 was *not run
   here* because IQ-TREE and MACSE are absent; `environment.yml` carries neither, there is no
   `docs/tools/` page for either, and no adapter exists. Under `docs/tools/README.md`'s role
   table, tree inference is `future` — a deferred branch with no module.

**What the tests below pin, so the proposal cannot quietly become "done" without the parts.**

These tests are `xfail(strict=True)`-free: they assert the *current* absence, so they pass today
and must be updated deliberately when the parts arrive. A test that asserted the module existed
could not be written before it does.

**The decision this proposal asks for.** The group's own route is PhyloSuite; CGAS's is
MAFFT → MACSE → IQ-TREE. These are different tools with different defaults, and a module must
pick one and say why. That is a scope decision for the maintainer, not something to infer from a
paper's citation list.
"""

from __future__ import annotations

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
MODULES = ROOT / "src" / "plastid" / "modules"
RESULT = ROOT / "src" / "plastid" / "result"


def test_the_phylogeny_module_exists_and_declares_its_tools() -> None:
    """M11.21 landed: the module exists, and it names the two tools it calls.

    The earlier version of this file asserted the module's *absence*. It is replaced rather than
    deleted so the history shows the change was deliberate: the module, its adapters, its
    inventory rows and this test arrived together.
    """
    from plastid.modules import registry

    spec = registry.get("phylogeny")
    assert spec.arity >= 2
    assert "mafft" in spec.citations and "iqtree" in spec.citations
    declared = {param.name for param in spec.params}
    assert {"iqtree", "model", "ultrafast_bootstrap", "seed"} <= declared


def test_the_supermatrix_is_written_as_an_artifact() -> None:
    """The alignment a tree step reads is a named artifact, not a private temporary.

    Before M11.21 no module exported an alignment: `diversity`/`popgen`/`kaks` built one in
    memory and exported tables only, so a tree module would have had to re-align — a second
    numeric path. `phylogeny` carries `supermatrix.fasta` as a Result artifact, so the payload
    writer puts it beside the tree and the matrix can be read back and aligned against.
    """
    text = (MODULES / "phylogeny.py").read_text(encoding="utf-8")
    assert '"supermatrix.fasta": matrix_text' in text
    assert '"tree.treefile": result.treefile' in text


def test_the_inference_tool_is_pinned_and_declared() -> None:
    """The third input, now present: a pinned tool with an adapter, a page and a role.

    The earlier test asserted no inference tool was declared. It is replaced by the positive
    check, because a tool in the environment without a declaration is what G-04 removed.
    """
    runtime = (ROOT / "environment.yml").read_text(encoding="utf-8")
    assert "iqtree=2.3.6" in runtime
    assert (ROOT / "docs" / "tools" / "iqtree.md").is_file()
    inventory = (ROOT / "docs" / "tools" / "README.md").read_text(encoding="utf-8")
    assert "external/iqtree.py" in inventory


def test_the_scope_page_says_what_the_tree_module_does_not_do() -> None:
    """The honest half: a tree module that does not choose an outgroup must say so.

    Asserting *"we can build a tree"* is not enough. The page has to keep saying that rooting is
    the user's choice and that a topology is not a verifiable number — otherwise a user reads
    "phylogeny" as "the toolkit decided my tree is right".
    """
    text = (ROOT / "docs" / "things_we_do_not_do.md").read_text(encoding="utf-8")
    assert "Choose an outgroup or root the tree for you" in text
    assert "not a verifiable number" in text


@pytest.mark.parametrize("part", ["alignment artifact", "supermatrix", "inference adapter"])
def test_the_proposal_lists_each_missing_part(part: str) -> None:
    """The proposal must keep naming all three, not shrink to whichever is convenient."""
    record = ROOT / ".dev" / "records" / "M11.21_tree_proposal_2026-10-28.md"
    if not record.is_file():
        pytest.skip("the M11.21 proposal record is not in this tree (.dev/ is maintainer-only)")
    text = record.read_text(encoding="utf-8").lower()
    assert part in text
