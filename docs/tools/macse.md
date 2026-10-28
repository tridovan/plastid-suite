# MACSE 2.07 — codon-aware CDS alignment (M11.21, opt-in)

**Why this file exists.** CLAUDE.md's invariant: a module that calls an external tool declares it
in the repository. `phylogeny` can align its per-gene CDS blocks with MAFFT (default) or with
MACSE (`--set aligner=macse`). MACSE aligns **nucleotide sequences through their amino-acid
translation**, so a frameshift is placed as a codon-level indel instead of three independent gaps.
This is the route CGAS's module 14 offers with `--macse`; here it is opt-in the same way, and the
default — the path every number in this repository was measured with — remains MAFFT.

**What it is.** MACSE (Multiple Alignment of Coding SEquences Accounting for Frameshifts and Stop
Codons), v2.07, a Java program by Ranwez et al. Cite:

> Ranwez V, Douzery EJP, Cambon C, Chantret N, Delsuc F. MACSE v2: Toolkit for the Alignment of
> Coding Sequences Accounting for Frameshifts and Stop Codons. *Mol Biol Evol* 35(10):2582–2584,
> 2018.

**Obtained, not bundled.** MACSE is under **CeCILL 2.1**, a free-software licence; the toolkit does
not redistribute the jar. The vendor URL recorded by CGAS's own INSTALL.md
(`https://bioweb.supagro.inra.fr/macse/releases/macse_v2.07.jar`) was measured on **28/10/2026**:
it redirects to `https://www.agap-ge2pop.org/macse/releases/macse_v2.07.jar` and answers **404**.
The pinned source is therefore the **bioconda** package:

```sh
# the package (noarch, 381 kB)
curl -fL -o macse-2.07.tar.bz2 \
  https://conda.anaconda.org/bioconda/noarch/macse-2.07-hdfd78af_0.tar.bz2
sha256sum macse-2.07.tar.bz2
# c68ac215cc94826c075c71ea8d9b00acc4bfed0248b95d49cc7edc68b59429e2

tar xjf macse-2.07.tar.bz2          # bin/macse, share/macse-2.07-0/macse_v2.07.jar
sha256sum share/macse-2.07-0/macse_v2.07.jar
# 96873a1465f1e1aa9d0c462d469eee6954d36e8639727f69f34a72fcb6583963

java -jar share/macse-2.07-0/macse_v2.07.jar -help
# This is MACSE V2.07 If you find MACSE useful, please cite: ...
```

`conda install -c bioconda macse=2.07` is the same package and puts a `macse` wrapper on `PATH`;
that line is in `environment.yml`. Either route works: the adapter takes the wrapper on `PATH`, or
`--set macse=/path/to/macse_v2.07.jar`.

**The release is checked, not assumed.** MACSE has no version flag; the only place it names its
release is the first line of `-help`, `This is MACSE V2.07`. The adapter runs `-help`, reads that
banner, and **refuses any release but 2.07** — the output alphabet and the alignment algorithm are
part of the result.

**The invocation the adapter uses**, one gene at a time:

```sh
java -jar macse_v2.07.jar -prog alignSequences \
  -seq input.fasta -out_NT nt.fasta -out_AA aa.fasta
```

**Two characters are not DNA.** MACSE writes a frameshift as `!` and a stop codon as `*` in the
nucleotide alignment. The adapter replaces both with `-` at the boundary, because the rest of the
toolkit reads ACGT (plus `N`/`-`) and a column holding `!` or `*` cannot be classified. The rule is
stated in the run's scope label ("`!` and `*` read as gaps") and in
`external/macse.py`; it is not silent.

**Measured on this machine, 28/10/2026** (`NC_008641.1`, `NC_016708.2`, `NC_028027.1`,
`NC_016736.1`; one gene, `rbcL`, 1,428–1,443 bp): **7.5 s**, exit 0, output alphabet `-ACGT`,
trailing gaps pad the shorter records. That is one gene; a full 4-genome run is **78 genes**, so
expect minutes, not seconds.

**Measured at 34 genomes, 28–29/10/2026** (34-record Iridaceae panel): with the default 600 s
per-gene limit the run stopped at gene 63 — 62 genes had taken ~63 s each and one exceeded the
limit (`TimeoutExpired`, manifest `failed`). With `--set macse_timeout=3600` the same run
**completed**: 78 CDS blocks, **70,131 sites**, MACSE **2 h 06**, IQ-TREE **4 min 52 s**, total
**2 h 11** on this machine, tree with mostly 100/100 support. The four-genome panel takes 22:55
for 70,110 sites, so the cost is the number of taxa (≈5.5×), not the number of sites. A
family-scale matrix needs `macse_timeout` above the 600 s default.

**Not vendored, deliberately:** a Java jar under CeCILL 2.1, obtained by the user — the same
reasoning as MUMmer, BLAST+ and IQ-TREE.
