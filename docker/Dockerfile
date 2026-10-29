# PlastidToolSuite — the release image (M9.2 · M9.3 · M9.4).
#
# **Who it is for.** Windows users: the toolkit does not run natively there, so this image *is*
# the Windows route (README → Install). It is also the reproducible environment for anyone who
# wants the pinned tools without conda on their own machine — the same versions CI installs, so a
# number measured here is the number measured there.
#
# **Versions are part of the result** (two MAFFT builds gave two Ts/Tv numbers
# from one genome). Every pin below is one the repository records — the runtime ones in
# `environment.yml`, Vmatch in `environment-check.yml` — and the build *prints* them: a version
# nobody can read is not a pin.
#
# **Licence (M9.3 — docs/licensing.md).** This image is built on the user's machine, from recorded
# URLs with recorded hashes, so what it contains is the user's own copies for their own academic
# work: MISA (no stated licence), Vmatch (academic-only, so not for commercial use) and table2asn
# (NCBI). **Do not publish a prebuilt image**: that would be redistribution, the one thing MISA's
# missing licence forbids — the project ships the recipe, not the image.
#
# **The golden numbers run at build time.** If they do not reproduce inside the image, the image
# does not build — the phase's own check (ROADMAP → Phase 9: one input, one set of numbers,
# three environments: WSL · CI · Docker).

FROM mambaorg/micromamba:2.0.5

# This image puts the base environment on PATH for `RUN` steps **only** when this ARG is declared.
# Without it the build dies at `RUN pip install …` with "pip: command not found" — measured
# 11/10/2026, on the first build of this file.
ARG MAMBA_DOCKERFILE_ACTIVATE=1

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates curl unzip \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# The pinned tool set — perl (MISA is a Perl script), MAFFT, MUMmer 4, PAML, BLAST+, ARAGORN,
# IQ-TREE (the `phylogeny` module's inference step, M11.21) and
# Vmatch (academic-only: the *second engine* the repeats count is checked against — optional in the
# toolkit, installed here so the checks that need it run rather than skip). Vmatch is a
# `cross-check` tool in `docs/tools/README.md`, not a runtime one: it is in this image because the
# tier-B golden compares the two engines, and a job that silently skipped would prove nothing.
# Nothing from `environment-future.yml` is here — a deferred branch has no business in a
# release image.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.12 \
        perl=5.32.1 \
        mafft=7.526 \
        mummer4=4.0.1 \
        paml=4.10.10 \
        blast=2.17.0 \
        aragorn=1.2.41 \
        iqtree=2.3.6 \
        vmatch=2.3.1 \
        r-base=4.3 \
        r-ggplot2 \
        r-dplyr \
        r-tidyr \
        r-reshape2 \
        r-rcolorbrewer \
        r-zoo \
        r-scales \
        r-pheatmap \
        r-ape \
    && micromamba clean -a -y

# **table2asn is deliberately NOT in this image.** It is the one tool that (a) is optional — the
# `.fsa` + `.tbl` route needs nothing (D3-04) — and (b) is a 56 MB fetch from NCBI that cannot come
# from conda. Measured 11/10/2026: through the Docker VM that download crawls (the host was pulling
# the same file at ~600 kB/s), and an optional tool must never be able to hang a user's build. So
# the image carries the tools the *analyses* need, and the `.sqn` step is what `docs/tools/
# table2asn.md` describes: fetch the binary on the host and mount it, or pass
# `--set table2asn=/path/to/table2asn`. Its single check is deselected by name below.

# MISA 2.1 is fetched, never baked in from our side: the toolkit stopped shipping it because its
# zip states no licence (D-q, docs/tools/misa.md). Same URL, same md5s — and the module re-reads
# the release from the script's own header, so a different MISA could not pass as this one.
ENV PTS_MISA_DIR=/home/mambauser/misa
RUN mkdir -p "$PTS_MISA_DIR" \
    && cd /tmp \
    && curl -sSL --retry 3 --max-time 300 -O \
         https://webblast.ipk-gatersleben.de/misa/misa_sourcecode_25082020.zip \
    && echo "cf9d3579f7caacd02f99830b329c7d9c  misa_sourcecode_25082020.zip" | md5sum -c - \
    && unzip -j misa_sourcecode_25082020.zip misa.pl -d "$PTS_MISA_DIR" \
    && echo "cfebc1d1a36fddfd6e0d897c96babd12  $PTS_MISA_DIR/misa.pl" | md5sum -c - \
    && rm misa_sourcecode_25082020.zip

WORKDIR /work
COPY --chown=$MAMBA_USER:$MAMBA_USER . /work
# `pytest` is named explicitly: pip does not read PEP 735 dependency-groups (dev). The optional
# extras are installed too (M11.8): without them the in-image golden skipped 61 tests — every
# figure, every .docx and every R-script check — and a build that skips half the suite proves
# nothing. They are optional for a *user*; they are not optional for this check.
RUN pip install --no-cache-dir '.[figures,docx,excel,stats]' pytest

# The pins, printed where they can be read.
RUN mafft --version \
    && nucmer --version | head -1 \
    && blastn -version | head -1 \
    && aragorn -h 2>&1 | grep -m1 "^ARAGORN" \
    && iqtree2 --version | head -1 \
    && perl -v | grep -E "^This is perl 5" \
    && python -c "import plastid; print('plastid', plastid.__version__)" \
    && ls -l "$PTS_MISA_DIR/misa.pl"

# The golden numbers, inside the container. `PTS_REQUIRE_TOOLS=1` turns a missing tool from a
# printed skip into a failure — the whole point of checking in here, and why Vmatch is installed
# above rather than left out. **One** check is deselected, **by name and with its reason printed**:
# table2asn's `.sqn` step, which this image does not carry (see the note above). A deselect that
# says what it dropped is not a silent skip; nothing else may skip.
# The pinned panels (`.panel/`) are not in the build context by design — they are gitignored,
# ~250 MB, and third-party data, so the panel half of the golden skips here and the number printed
# is the vendored-fixture half only. To run the whole suite, including the panel:
#
#   docker run --rm -v "$PWD:/work" -w /work pts-local:golden \
#     sh -c 'PTS_REQUIRE_TOOLS=1 python -m pytest tests/golden -q -rs \
#       --deselect tests/golden/test_tier_b.py::test_the_sqn_step_runs_and_ncbi_reads_the_table'
#
# Measured 28/10/2026: 55 skips inside the build, every one of them "`.panel/` is not in this
# tree"; **0 skips** with the repository mounted. And the skip count is not silent: the suite
# prints every reason, which is why this note exists.
RUN echo "table2asn is not in this image: its one optional .sqn check is deselected by name" \
    && echo "panel checks skip here: .panel/ is mounted at run time, not built in" \
    && PTS_REQUIRE_TOOLS=1 python -m pytest tests/golden -q -rs \
         --deselect tests/golden/test_tier_b.py::test_the_sqn_step_runs_and_ncbi_reads_the_table

CMD ["plastid", "--help"]
