UV ?= uv

.PHONY: env test golden test-all lint fmt check

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
