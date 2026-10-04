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
        {"annotator": "pga", "reference": str(reference), "pga": str(_stub_pga(tmp_path))},
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
    # Ten features carry a second name: ycf3/ycf4/psbN/clpP and ndhL (the synonym this record
    # uses for ndhI, 2026-10-03); atpI, psaI and psbB (the Ka/Ks run showed NC_024936.1 spells
    # those three loci atpL, psaL and bspB); and both IR copies of trnI-GAU (that record spells
    # the locus trnl-GUA, and the sequences are byte-identical). Every entry's evidence is in
    # data/gene_aliases.yaml.
    assert summary["n_aliased_names"] == 10

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
        annotation.run(
            read_genbank(data_dir / "NC_050833.1.gb"), {"annotator": "pga"}, _ctx(tmp_path)
        )


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _tiny_genbank(*, with_genes: bool, partial: bool = False) -> str:
    """A minimal record: one 600 bp sequence, optionally with two calls of one locus."""
    sequence = ("atg" + "ttc" * 190 + "taa" + "a" * 24)[:600]  # 600 bp, CDS 1..576
    lines = [
        "LOCUS       SYN_ANNOT                600 bp    DNA     linear   PLN 03-OCT-2026",
        "FEATURES             Location/Qualifiers",
        "     source          1..600",
    ]
    if with_genes:
        location = "<1..576" if partial else "1..576"
        for name in ("clpP", "clpP1"):
            lines += [
                f"     gene            {location}",
                f'                     /gene="{name}"',
                f"     CDS             {location}",
                f'                     /gene="{name}"',
                '                     /product="clp protease proteolytic subunit"',
                "                     /transl_table=11",
            ]
    origin = [
        f"{start:>9} "
        + " ".join(
            sequence[i : i + 10] for i in range(start - 1, min(start - 1 + 60, len(sequence)), 10)
        )
        for start in range(1, len(sequence) + 1, 60)
    ]
    lines += ["ORIGIN", *origin, "//"]
    return "\n".join(lines) + "\n"


def test_two_calls_of_one_locus_are_reported_not_merged(tmp_path: Path, monkeypatch) -> None:
    """Measured with PGA and seven references (2026-10-03): one clpP locus came back twice, as
    clpP and clpP1, because some references use the classic names and some the functional ones.
    The two rows must both exist and the duplication must be named — merging would hide a real
    property of the annotation, and dropping one would hide the other name."""
    from plastid.modules import annotation

    source = _write(tmp_path / "annotated.gb", _tiny_genbank(with_genes=True))
    monkeypatch.setenv("STUB_GENBANK", str(source))
    reference = tmp_path / "refs"
    reference.mkdir()
    _write(reference / "ref.gb", _tiny_genbank(with_genes=True))

    result = annotation.run(
        read_genbank(source),
        {"annotator": "pga", "reference": str(reference), "pga": str(_stub_pga(tmp_path))},
        _ctx(tmp_path),
    )

    assert len(result.tables["gene_set"].rows) == 2  # both calls kept
    duplicates = result.tables["duplicate_loci"].rows
    assert [row["canonical"] for row in duplicates] == ["clpP1"]
    assert duplicates[0]["names"] == "clpP, clpP1"
    assert any("annotated more than once" in note for note in result.notes)


def test_the_two_inverted_repeat_copies_are_not_reported_as_duplicates(
    tmp_path: Path, monkeypatch
) -> None:
    """A plastome carries two copies of every IR gene, under the same name. Counting those as
    duplicates buried the finding that matters: measured 2026-10-03, 25 "duplicates" of which 21
    were the two copies of an IR gene."""
    from plastid.modules import annotation

    text = _tiny_genbank(with_genes=False).replace(
        "     source          1..600",
        "     source          1..600\n"
        '     gene            1..576\n                     /gene="ndhB"\n'
        '     CDS             1..576\n                     /gene="ndhB"\n'
        '                     /product="NADH dehydrogenase subunit B"\n'
        "                     /transl_table=11\n"
        '     gene            1..576\n                     /gene="ndhB"\n'
        '     CDS             1..576\n                     /gene="ndhB"\n'
        '                     /product="NADH dehydrogenase subunit B"\n'
        "                     /transl_table=11",
    )
    source = _write(tmp_path / "ircopies.gb", text)
    monkeypatch.setenv("STUB_GENBANK", str(source))
    reference = tmp_path / "refs_ir"
    reference.mkdir()
    _write(reference / "ref.gb", _tiny_genbank(with_genes=True))

    result = annotation.run(
        read_genbank(source),
        {"annotator": "pga", "reference": str(reference), "pga": str(_stub_pga(tmp_path))},
        _ctx(tmp_path),
    )
    assert len(result.tables["gene_set"].rows) == 2
    assert "duplicate_loci" not in result.tables


def test_a_reference_with_partial_coordinates_is_refused_by_name(tmp_path: Path) -> None:
    """Measured 2026-10-03: PGA cannot parse `<`/`>` — it warns, annotates nothing, and exits 0.
    A silent empty annotation is the failure this check turns into a named one."""
    from plastid.external.pga import run_pga

    reference = tmp_path / "refs_partial"
    reference.mkdir()
    _write(reference / "partial.gb", _tiny_genbank(with_genes=True, partial=True))

    with pytest.raises(ToolError, match="partial coordinates"):
        run_pga(
            read_genbank(_write(tmp_path / "target.gb", _tiny_genbank(with_genes=False))),
            _ctx(tmp_path),
            reference_dir=reference,
            seq_id="SYN_ANNOT",
            script=_stub_pga(tmp_path),
        )


def test_an_annotation_with_no_genes_is_a_named_error_not_an_empty_report(
    tmp_path: Path, monkeypatch
) -> None:
    """PGA reports success while annotating nothing (measured). A run that produced no gene is a
    failed annotation, not a result with zeros."""
    from plastid.external.pga import run_pga

    empty = _write(tmp_path / "empty.gb", _tiny_genbank(with_genes=False))
    monkeypatch.setenv("STUB_GENBANK", str(empty))
    reference = tmp_path / "refs_ok"
    reference.mkdir()
    _write(reference / "ref.gb", _tiny_genbank(with_genes=True))

    with pytest.raises(ToolError, match="no gene features"):
        run_pga(
            read_genbank(empty),
            _ctx(tmp_path),
            reference_dir=reference,
            seq_id="SYN_ANNOT",
            script=_stub_pga(tmp_path),
        )


_FLAKY_STUB = r"""#!/usr/bin/env perl
use strict;
use warnings;
my ($out, $target) = ('', '');
my @args = @ARGV;
while (@args) {
    my $flag = shift @args;
    $out = shift @args if $flag eq '-o';
    $target = shift @args if $flag eq '-t';
}
mkdir $out;
# the caller owns the counter: the stub only reads it, so there is no doubt about who created it
open(my $cfh, '<', $ENV{STUB_COUNTER}) or die "no counter\n";
my $attempt = <$cfh>;
close $cfh;
$attempt++;
open(my $cfw, '>', $ENV{STUB_COUNTER}) or die $!;
print $cfw $attempt;
close $cfw;
opendir(my $dh, $target) or die $!;
my @fasta = grep { /\.fasta$/ } readdir($dh);
closedir($dh);
(my $stem = $fasta[0]) =~ s/\.fasta$//;
my $source = ($attempt == 1) ? $ENV{STUB_EMPTY} : $ENV{STUB_GENBANK};
open(my $in, '<', $source) or die $!;
my $text = do { local $/; <$in> };
close $in;
open(my $outfh, '>', "$out/$stem.gb") or die $!;
print $outfh $text;
close $outfh;
"""


def _flaky_pga(tmp_path: Path) -> Path:
    script = tmp_path / "PGA_flaky.pl"
    script.write_text(_FLAKY_STUB, encoding="utf-8")
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return script


def _prepare(tmp_path: Path, monkeypatch) -> tuple[Path, Path]:
    annotated = _write(tmp_path / "genes.gb", _tiny_genbank(with_genes=True))
    empty = _write(tmp_path / "none.gb", _tiny_genbank(with_genes=False))
    counter = _write(tmp_path / "attempts.count", "0")
    monkeypatch.setenv("STUB_GENBANK", str(annotated))
    monkeypatch.setenv("STUB_EMPTY", str(empty))
    monkeypatch.setenv("STUB_COUNTER", str(counter))
    reference = tmp_path / "refs_flaky"
    reference.mkdir()
    _write(reference / "ref.gb", _tiny_genbank(with_genes=True))
    return empty, reference


def test_a_first_attempt_that_finds_nothing_is_retried_and_the_report_says_so(
    tmp_path: Path, monkeypatch
) -> None:
    """Measured 2026-10-03: PGA's temporary names can contain '$', which the shell expands away
    on its way to blastn — about one run in five produces nothing and still exits 0. One retry
    turns that into a working run, and the report admits the retry happened."""
    from plastid.external.pga import run_pga

    empty, reference = _prepare(tmp_path, monkeypatch)
    result = run_pga(
        read_genbank(empty),
        _ctx(tmp_path),
        reference_dir=reference,
        seq_id="SYN_ANNOT",
        script=_flaky_pga(tmp_path),
    )
    print("COUNTER after:", (tmp_path / "attempts.count").read_text())
    assert result.attempts == 2
    assert [gene.name for gene in result.genome.genes if gene.kind == "CDS"]


def test_two_empty_attempts_are_a_named_error(tmp_path: Path, monkeypatch) -> None:
    """A stub that always writes the empty record: the retry must not hide a real failure."""
    from plastid.external.pga import run_pga

    empty, reference = _prepare(tmp_path, monkeypatch)
    monkeypatch.setenv("STUB_GENBANK", str(empty))  # both attempts come back empty
    with pytest.raises(ToolError, match="two attempts"):
        run_pga(
            read_genbank(empty),
            _ctx(tmp_path),
            reference_dir=reference,
            seq_id="SYN_ANNOT",
            script=_flaky_pga(tmp_path),
        )
