UV ?= uv

.PHONY: env test golden test-all lint fmt check corpus-bar

env:
	$(UV) sync

test:
	$(UV) run pytest -q -m "not tier_b"

golden:
	$(UV) run pytest -q tests/golden -rs

test-all:
	$(UV) run pytest -q -rs

lint:
	$(UV) run ruff check .
	$(UV) run ruff format --check .

fmt:
	$(UV) run ruff format .

check: lint test golden

# P3B-10's bar, measured: 34 leave-one-out runs of the pinned corpus, each classified against the
# published record and against PGA's reading. Needs `.panel/` (gitignored), BLAST+ and ARAGORN —
# `make corpus-bar TOOLS=/path/to/bin` when they are not on PATH. Minutes, not seconds; not part
# of `check`.
corpus-bar:
	$(UV) run python tools/corpus_bar.py all $(if $(TOOLS),--tools $(TOOLS),)
