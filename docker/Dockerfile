# Skeleton — the release image is built and verified in Phase 9 (ROADMAP.md M9.2, M9.4).
# Windows users run the toolkit through this image; macOS stays experimental.
FROM python:3.12-slim

# MISA 2.1 (vendored, md5-pinned in src/plastid/vendor/misa/SOURCE.json) needs perl.
RUN apt-get update \
    && apt-get install -y --no-install-recommends perl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /work
COPY . /work
# pytest is installed explicitly: pip does not read PEP 735 dependency-groups (dev).
RUN pip install --no-cache-dir -e . pytest

# Smoke check: the golden numbers must reproduce inside the container.
RUN python -m pytest tests/golden -q -rs

CMD ["plastid", "--version"]
