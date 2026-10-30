UV ?= uv

.PHONY: env test golden acceptance test-all lint fmt check corpus-bar

env:
	$(UV) sync

test:
	$(UV) run pytest -q -m "not tier_b and not panel"

# The golden numbers on the golden pair and the pinned-tool checks — the set CI runs on every push.
golden:
	$(UV) run pytest -q tests/golden -rs -m "not panel" --durations=25

# Every test that reads `.panel/` — 80 RefSeq plastomes in 10 families, 30 PGA re-annotations,
# rotated FASTA, the pinned reference sets. `.panel/` is gitignored, so CI never runs this tier;
# run it when closing a milestone. A skip is not acceptance, so a missing panel or tool fails here.
acceptance:
	@test -d .panel || { echo "make acceptance: .panel/ is not in this tree"; exit 1; }
	PTS_REQUIRE_TOOLS=1 PTS_REQUIRE_PANEL=1 $(UV) run pytest -q -rs -m panel --durations=25

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
