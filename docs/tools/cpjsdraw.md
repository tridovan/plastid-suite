# CPJSdraw — declared here, obtained separately (the source states no licence)

**Why this file exists.** The toolkit never bundles a tool it can call, and this one must not be an
exception: the CPJSdraw tree carries **no licence file** and its `readme.txt` states none. Measured
17/10/2026: no `LICENSE`, no `COPYING` in the tree, and GitHub's API answers `license: null` for
`xul962464/CPJSdraw`. Under the rule that made MISA a declaration instead of a vendored script
(P9-01) — *no licence means all rights reserved, not "free to take"* — CPJSdraw is **cloned by the
user**, never copied into this repository or the wheel.

> The integration bundle this adapter ports contradicts itself here: its `README_tich_hop.md` says
> copying the tree into `vendor/` is "valid because CPJSdraw is MIT", while its own
> `docs_tools_mau/cpjsdraw.md` and `PATCH_DA_AP_DUNG.md` both record that the tree has no licence.
> The tree is what was checked; there is no MIT.

**What it is.** A standalone Perl program that draws chloroplast **junction sites** (JLB/JSB/JSA/JLA)
for one or many GenBank records: each record's LSC/IRb/SSC/IRa is supplied (or found with MUMmer4),
and the flanking genes are drawn to scale on both sides of every junction. It is the tool the
Iridaceae manuscript's Figure 2 was drawn with.

**Not CPStools.** CPStools (Xwb7533) is a different package by a different author. Do not conflate
them, in code or in a manuscript's methods.

**Cite it — the citable paper is PeerJ, not the GitHub page:**

> Li H, Guo Q, Xu L, Gao H, Liu L, Zhou X. **CPJSdraw: analysis and visualization of junction sites
> of chloroplast genomes.** *PeerJ* 11:e15326 (2023). doi:10.7717/peerj.15326

*(The repository's `readme.txt` gives exactly this citation. A manuscript citing "Xu, W. (2024),
[Software], GitHub" has the wrong year, the wrong first author and a non-archival source — the
bundle's own methods section did, and it is corrected before it reaches a journal.)*

**Source, pinned (measured 17/10/2026):**

| | |
|---|---|
| URL | `https://github.com/xul962464/CPJSdraw.git` |
| The script's own version | `Version: 0.0.1`, read from `CPJSdraw.pl`'s header (it has no version flag) |
| `bin/CPJSdraw.pl` **as patched here** | 52,947 B · md5 `ba41be6002bfd0b797907d3290b34419` |
| Dependencies | Perl + the `SVG` module (`apt install libsvg-perl`); **MUMmer4** (`nucmer`, `show-coords` on `PATH`) only when a record's region string is not supplied — this toolkit always supplies it |

**Install line (verified 17/10/2026):**

```sh
git clone https://github.com/xul962464/CPJSdraw.git
chmod a+x CPJSdraw/script/svg_kit/buildInFont      # readme.txt asks for this
sudo apt-get install -y perl libsvg-perl librsvg2-bin
```

**It wants `gene` features, and this toolkit writes none.** Measured 17/10/2026: handed a record
without one, CPJSdraw stops with `has no gene info(no primary gene tag)` and writes no SVG. A
GenBank record from NCBI carries both `gene` and `CDS` features and works; a record **this toolkit
annotated** carries `CDS`/`tRNA`/`rRNA` only, so `plastid annotate` output cannot be drawn as it
stands. The refusal is loud — an exit-3 with the tool's own message — rather than a figure that is
quietly missing, which is why it is written down here instead of worked around.

**Nothing is probed with a flag.** `CPJSdraw.pl` is a Perl script: passing it an unknown flag makes
it run itself with default arguments. Its version is read out of its own header, and the
interpreter's version is asked of Perl (`perl -e 'print $^V'`).

**SVG → PNG uses `rsvg-convert`, not the bundled batik.** CPJSdraw calls `script/svg_kit/svg2xxx`,
which needs Java and the `batik-1.7` it ships; on Java 17+ that fails, because batik 1.7 uses the
removed Security Manager. Measured in the bundle: `rsvg-convert -w 3200 <svg> -o <png>` produces the
same picture, and it is what the manuscript's PNGs were made with.

## Two patches — the figure cannot be reproduced from a pristine clone

`src/plastid/external/cpjsdraw.py` **refuses** a script whose md5 is not the one above, naming these
hunks. The refusal is deliberate: the second patch is a correctness fix, and without it the figure
is wrong in a way nobody notices.

**1. A bolder face** — `PPI` 300 → 600, larger fonts, thicker strokes, `font-weight='bold'` on every
text, arrow heads ±3 → ±5. Check: `grep -c bold bin/CPJSdraw.pl` → **17** on the patched file,
**0** on a clone.

**2. A compound IRa** — the fix that matters. The regex reads `IRa:(\d+)-(\d+)` only, so when the
genome origin falls *inside* IRa (the copy is written `join(1..4,126609..152777)`) the part across
the origin is dropped and IRa is drawn short. Two edits:

```perl
# around line 136, and again in sub format_js_info (around line 1081):
- /LSC:(\d+)-(\d+),IRb:(\d+)-(\d+),SSC:(\d+)-(\d+),IRa:(\d+)-(\d+)/i
+ /LSC:(\d+)-(\d+),IRb:(\d+)-(\d+),SSC:(\d+)-(\d+),IRa:(\d+)-(\d+)(?:;(\d+)-(\d+))?/i

# and immediately after the match in format_js_info, BEFORE the for my $i(1..$#pos) loop:
+ if(defined $pos[8] && defined $pos[9]){
+     $pos[7] = $pos[9] + $genome_len;   # IRa = its main run plus the part across the origin
+     splice(@pos, 8, 2);                # drop the two wrap coordinates
+ }
```

Shortfall without it, measured on the five records this project drew that wrap: **Chasmanthe
floribunda 1 bp · Trimezia steyermarkii 4 bp · Isophysis tasmanica 7 bp · Freesia laxa 36 bp ·
Olsynium douglasii 42 bp**.

*(The bundle's methods section names *Iris hermona* among the compound-IRa species. The genome set
it drew Figure 2 from has **Isophysis tasmanica** in that place, wrapping 7 bp, and Iris hermona's
copies do not wrap. Correct the list before the manuscript goes anywhere else.)*
