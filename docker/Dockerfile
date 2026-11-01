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
        macse=2.07 \
        openjdk=17 \
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

# The assembly tools and table2asn each take an environment of their own: GetOrganelle pins a BLAST
# older than the analyses', and table2asn cannot be solved beside R (docs/tools/table2asn.md). They
# come last on PATH, so the analyses keep their pinned BLAST. This layer sits above the COPY so a
# source change does not pay for it again. GetOrganelle's and Oatk's databases are fetched here, from
# the releases docs/tools/getorganelle.md and docs/tools/oatk.md record: someone on Windows has no
# prompt to fetch them from, and 140 MB is little beside the 3 GB of tools that need them.
COPY environment-assembly.yml /tmp/environment-assembly.yml
RUN micromamba create -y -n assembly -f /tmp/environment-assembly.yml \
    && micromamba create -y -n table2asn -c conda-forge -c bioconda table2asn=1.29.324 \
    && micromamba clean -a -y
ENV PATH=$PATH:/opt/conda/envs/assembly/bin:/opt/conda/envs/table2asn/bin
ENV GETORG_PATH=/home/mambauser/db/getorganelle
ENV OATK_DB=/home/mambauser/db/oatk/v20230921
# The assembly environment first, so this database is built by the BLAST that reads it. The archive
# holds every clade, 210 MB; the plant profiles are the 61 MB a plastome needs.
RUN mkdir -p "$GETORG_PATH" "$OATK_DB" \
    && PATH=/opt/conda/envs/assembly/bin:$PATH \
       get_organelle_config.py -a embplant_pt,embplant_mt --config-dir "$GETORG_PATH" \
    && curl -sSL --retry 3 --max-time 300 -o /tmp/oatkdb.tar.gz \
         https://github.com/c-zhou/OatkDB/archive/75e8db0ac4a7d508a9a518d900876003ceb70737.tar.gz \
    && tar -xzf /tmp/oatkdb.tar.gz -C "$OATK_DB" --strip-components=2 \
         --wildcards '*/v20230921/embryophyta_*' \
    && rm /tmp/oatkdb.tar.gz \
    && ls "$OATK_DB"

# What CPJSdraw needs to run. CPJSdraw itself is not here: its source states no licence and the
# figure needs a patched build (docs/tools/cpjsdraw.md), so the user names their own copy in the
# mounted folder, on the Export page or with --cpjsdraw.
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends libsvg-perl librsvg2-bin \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

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
    && ls -l "$PTS_MISA_DIR/misa.pl" \
    && get_organelle_from_reads.py --version \
    && oatk --version \
    && table2asn -version \
    && plastid tools

# The golden numbers, inside the container. `PTS_REQUIRE_TOOLS=1` turns a missing tool from a
# printed skip into a failure, the whole point of checking in here.
# The pinned panels (`.panel/`) are not in the build context by design: they are gitignored,
# ~250 MB, and third-party data, so the panel half of the golden skips here and the number printed
# is the vendored-fixture half only. To run the whole suite, including the panel:
#
#   docker run --rm -v "$PWD:/work" -w /work pts-local:golden \
#     sh -c 'PTS_REQUIRE_TOOLS=1 python -m pytest tests/golden -q -rs'
RUN echo "panel checks skip here: .panel/ is mounted at run time, not built in" \
    && PTS_REQUIRE_TOOLS=1 python -m pytest tests/golden -q -rs

# The container listens on every address it has, because a published port arrives from outside it.
# `docker/plastid-ui.bat` publishes that port to 127.0.0.1 on the host alone, and the session token,
# the Host check and the Origin check stand exactly as they do outside a container.
ENV PLASTID_IN_CONTAINER=1
EXPOSE 26214
CMD ["plastid", "ui", "--workspace", "/work", "--host", "0.0.0.0", "--no-browser"]
