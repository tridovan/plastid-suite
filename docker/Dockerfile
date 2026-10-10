# PlastidToolSuite — the release image (M9.2 · M9.3 · M9.4).
#
# **Who it is for.** Windows users: the toolkit does not run natively there, so this image *is*
# the Windows route (README → Install). It is also the reproducible environment for anyone who
# wants the pinned tools without conda on their own machine — the same versions CI installs, so a
# number measured here is the number measured there.
#
# **Versions are part of the result** (PROJECT.md §9.B-14: two MAFFT builds gave two Ts/Tv numbers
# from one genome). Every pin below is the one `environment.yml` and §7.8 record, and the build
# *prints* them: a version nobody can read is not a pin.
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

# The pinned tool set — perl (MISA is a Perl script), MAFFT, MUMmer 4, PAML, BLAST+, ARAGORN, and
# Vmatch (academic-only: the *second engine* the repeats count is checked against — optional in the
# toolkit, installed here so the checks that need it run rather than skip).
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.12 \
        perl=5.32.1 \
        mafft=7.526 \
        mummer4=4.0.1 \
        paml=4.10.10 \
        blast=2.17.0 \
        aragorn=1.2.41 \
        vmatch=2.3.1 \
    && micromamba clean -a -y

# table2asn cannot come from conda (docs/tools/table2asn.md), so it is fetched here — **best
# effort**, because it is the one *optional* tool in the image (the `.fsa` + `.tbl` route needs
# nothing) and the one download a locked-down network may not reach. Measured 11/10/2026: through
# the Docker VM this download can crawl, and an optional tool must never be able to block a user's
# build. A **corrupted** download is still fatal — the md5 check aborts the build.
USER root
RUN set -e; \
    if curl -sSL --retry 3 --max-time 300 -o /tmp/table2asn.gz \
         https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/linux64.table2asn.gz; then \
      gunzip -f /tmp/table2asn.gz; \
      echo "8b5b6cbede9d38b47c7c4cae29a25a0a  /tmp/table2asn" | md5sum -c -; \
      install -m 0755 /tmp/table2asn /usr/local/bin/table2asn; \
      rm -f /tmp/table2asn; \
    else \
      echo "table2asn not fetched (network) — only the optional .sqn step needs it"; \
    fi
USER $MAMBA_USER

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
# `pytest` is named explicitly: pip does not read PEP 735 dependency-groups (dev).
RUN pip install --no-cache-dir . pytest

# The pins, printed where they can be read.
RUN mafft --version \
    && nucmer --version | head -1 \
    && blastn -version | head -1 \
    && aragorn -h 2>&1 | grep -m1 "^ARAGORN" \
    && perl -v | grep -E "^This is perl 5" \
    && python -c "import plastid; print('plastid', plastid.__version__)" \
    && ls -l "$PTS_MISA_DIR/misa.pl"

# The golden numbers, inside the container. `PTS_REQUIRE_TOOLS=1` turns a missing tool from a
# printed skip into a failure — the whole point of checking in here. The single exception is
# table2asn, and it is **deselected by name, with the reason printed**: it is optional (D3-04), it
# produces no §6 golden number, and it is the one tool this image may not have fetched.
RUN if command -v table2asn >/dev/null 2>&1; then \
      echo "table2asn present — every check runs, and a missing tool is a failure"; \
      PTS_REQUIRE_TOOLS=1 python -m pytest tests/golden -q -rs; \
    else \
      echo "table2asn absent — its one optional check is deselected by name; nothing else may skip"; \
      PTS_REQUIRE_TOOLS=1 python -m pytest tests/golden -q -rs \
        --deselect tests/golden/test_tier_b.py::test_the_sqn_step_runs_and_ncbi_reads_the_table; \
    fi

CMD ["plastid", "--help"]
