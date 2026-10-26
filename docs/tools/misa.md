# MISA 2.1 — declared here, and no longer bundled (D-q, 11/10/2026)

**Why this file exists.** CLAUDE.md: the moment a module uses an external tool, that tool is
declared in the repo — in `environment.yml`, in `vendor/` beside a `SOURCE.json`, or **here**. MISA
sat in the second slot from Phase 2 until 11/10/2026. It cannot stay there: **the source zip states
no licence** — it ships no licence file, and the IPK page only asks for the citation — and
redistributing a work is the one act a missing licence forbids (no licence means all rights
reserved, not "free to take"). So the script is **obtained, not shipped**, and this file is its
third-slot declaration. The rule, and the inventory it belongs to: `docs/licensing.md`.

**What it is.** MISA (MIcroSAtellite identification tool) 2.1 — Thomas Thiel & Sebastian Beier, IPK
Gatersleben. It is the SSR caller `ssr` counts with. **Cite it:** Beier S, Thiel T, Munch T,
Scholz U, Mascher M. *MISA-web: a web server for microsatellite prediction.* Bioinformatics
33:2583–2585 (2017). doi:10.1093/bioinformatics/btx198.

**Source, pinned** (recorded 02/10/2026 from the copy that was vendored; kept after it left the
tree, because a pin outlives the file):

| | |
|---|---|
| URL | `https://webblast.ipk-gatersleben.de/misa/misa_sourcecode_25082020.zip` |
| Zip | 4,946 B · md5 `cf9d3579f7caacd02f99830b329c7d9c` |
| `misa.pl` | 20,148 B · md5 `cfebc1d1a36fddfd6e0d897c96babd12` · sha256 `4712f29a0c57ff4aee1cc248cba20bcc879ef5d3b2ea6f5f1162cbc391c59abb` |
| Release | the script's own header: `###Release date: 25/08/20 (version 2.1)` |

**Install line** (verified 11/10/2026 on the dev WSL machine — the copy the toolkit used to ship
hashes to exactly this):

```sh
curl -sSLO https://webblast.ipk-gatersleben.de/misa/misa_sourcecode_25082020.zip
md5sum misa_sourcecode_25082020.zip      # cf9d3579f7caacd02f99830b329c7d9c
unzip -j misa_sourcecode_25082020.zip misa.pl -d "$HOME/misa"
md5sum "$HOME/misa/misa.pl"              # cfebc1d1a36fddfd6e0d897c96babd12
```

Then either put `misa.pl` on `PATH`, or point the module at it:

```sh
plastid run --module ssr --input <genome>.gb --out .runs --dataset <name> \
  --set misa="$HOME/misa/misa.pl"
```

**The release is checked, not assumed.** MISA has no version flag; `external/misa.py` reads
`(version 2.1)` out of the script's first lines and **refuses any other release by name** — the
same rule `repeat_match_version` applies to MUMmer, and for the same reason: the golden SSR
numbers (**37 = mono 31 · di 4 · tri 2**) were measured with 2.1, and a different caller's counts
are a different number.

**`misa.ini` is the toolkit's own file** (`src/plastid/data/misa.ini`), not the zip's: it writes
down MISA's **default** definitions (unit sizes `1-10 2-6 3-5 4-5 5-5 6-5`, interruptions `100`,
`GFF: true`). Those are the tool's documented defaults and the project uses exactly them —
re-measured 10/10/2026: at the defaults our table reproduces the published answer key,
while the thresholds that paper *states* in its Methods do not.

**Not vendored, deliberately** — the licence paragraph above. The consequence, stated so it is not
a surprise: `make golden`'s SSR test **skips with a named reason** on a machine without `misa.pl`
(and fails under `PTS_REQUIRE_TOOLS=1`), exactly like MAFFT, BLAST+ and ARAGORN. The release's
Docker image obtains it at build time from the recorded URL with the md5 checked (M9.2) — a
locally-built image is the user's own copy; **publishing a prebuilt image would be redistribution**
and needs IPK's written permission, which is the open ask on the board.

**Measured while unbundling** (11/10/2026): the wheel `uv build --wheel` produced before this
change carried **146 files**, `plastid/vendor/misa/misa.pl` among them — a bundled third-party work
whose licence is not stated. That measurement is what D-q answered.
