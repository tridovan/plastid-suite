"""PGA annotation (P3-03): the wiring around the annotator, and what it refuses to guess.

PGA is not installed here, so these tests exercise the two halves that do not need it: the named
error when the tool is missing, and the path from "the annotator wrote a GenBank file" to "the
gene set went through the same checks as a downloaded record". The second half runs against a
**stub** in place of PGA — it tests our code, not PGA's annotation quality, and it says so.
"""

from __future__ import annotations

import stat
from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank
from plastid.external.pga import run_pga
from plastid.run.context import RunContext
from plastid.run.tools import ToolError, which

_STUB = """#!/usr/bin/env perl
use strict;
use warnings;
my ($out, $target) = ('', '');
my @args = @ARGV;
while (@args) {
    my $flag = shift @args;
    $out = shift @args if $flag eq '-o';
    $target = shift @args if $flag eq '-t';
}
die "stub: no -o\\n" unless $out;
die "stub: no -t\\n" unless $target;
mkdir $out;
my $source = $ENV{STUB_GENBANK} or die "stub: STUB_GENBANK not set\\n";
opendir(my $dh, $target) or die $!;
my @fasta = grep { /\\.fasta$/ } readdir($dh);
closedir($dh);
die "stub: no target fasta\\n" unless @fasta;
(my $stem = $fasta[0]) =~ s/\\.fasta$//;
open(my $in, '<', $source) or die $!;
my $text = do { local $/; <$in> };
close $in;
open(my $outfh, '>', "$out/$stem.gb") or die $!;
print $outfh $text;
close $outfh;
"""


def _ctx(tmp_path: Path) -> RunContext:
    from plastid.run.layout import RunLayout

    layout = RunLayout.create(tmp_path / "runs", "annotation-test", fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id="annotation", module_order=1)


def _stub_pga(tmp_path: Path) -> Path:
    script = tmp_path / "PGA.pl"
    script.write_text(_STUB, encoding="utf-8")
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return script


def test_a_missing_pga_is_a_named_error_that_says_how_to_get_it(
    data_dir: Path, tmp_path: Path
) -> None:
    if which("PGA.pl") is not None:
        pytest.skip("PGA.pl is on PATH here — this test is about the machine that has none")
    genome = read_genbank(data_dir / "NC_050833.1.gb")
    with pytest.raises(ToolError, match="PGA.pl not found"):
        run_pga(genome, _ctx(tmp_path), reference_dir=tmp_path, seq_id="NC_050833.1")


def test_a_window_into_a_file_that_is_not_pga_pl_is_refused(data_dir: Path, tmp_path: Path) -> None:
    genome = read_genbank(data_dir / "NC_050833.1.gb")
    with pytest.raises(ToolError, match="pass --set pga="):
        run_pga(
            genome,
            _ctx(tmp_path),
            reference_dir=tmp_path,
            seq_id="NC_050833.1",
            script=tmp_path / "nope.pl",
        )


def test_the_module_analyses_whatever_the_annotator_wrote(
    data_dir: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The whole path with a stub in PGA's place: our argv, our read-back, our checks. The stub
    copies the golden genome's own GenBank, so the numbers below are the golden's — which is
    exactly what a working annotator is expected to reproduce (P3-03's Done-when)."""
    from plastid.modules import annotation

    golden_path = data_dir / "NC_050833.1.gb"
    monkeypatch.setenv("STUB_GENBANK", str(golden_path))
    reference = tmp_path / "refs"
    reference.mkdir()
    (reference / golden_path.name).write_bytes(golden_path.read_bytes())

    result = annotation.run(
        read_genbank(golden_path),
        {"reference": str(reference), "pga": str(_stub_pga(tmp_path))},
        _ctx(tmp_path),
    )

    summary = result.tables["annotation_summary"].rows[0]
    assert (summary["n_genes"], summary["n_cds"], summary["n_trna"], summary["n_rrna"]) == (
        132,
        86,
        38,
        8,
    )
    assert summary["n_orf_failed"] == 0
    assert summary["n_unclassified"] == 0
    assert summary["n_aliased_names"] == 4  # ycf3/ycf4/psbN/clpP carry a second name

    status = {row["status"]: row["n"] for row in result.tables["orf_status"].rows}
    assert status["ok"] == 81 and status["alternative_start"] == 5

    # the raw annotator output travels with the result, so a number can be traced back to it
    assert "annotation.gb" in result.artifacts
    assert "NC_050833.1" in result.artifacts["annotation.gb"]

    # the gene set is grouped and printed with both names, like every other table in a run
    genes = {row["gene"] for row in result.tables["gene_set"].rows}
    assert "pafI (ycf3)" in genes
    assert len(result.tables["gene_set"].rows) == 132


def test_the_module_names_the_missing_reference_directory(data_dir: Path, tmp_path: Path) -> None:
    from plastid.modules import annotation

    with pytest.raises(ValueError, match="reference="):
        annotation.run(read_genbank(data_dir / "NC_050833.1.gb"), {}, _ctx(tmp_path))
