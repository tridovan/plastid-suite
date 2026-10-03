# Packaged data

Resource files shipped inside the package (not a Python package — no `__init__.py`).

- `defaults.yaml` — central default configuration (arrives in P1-05)
- `citations.yaml` — citation registry: tool → reference (arrives in P1-07)

Access from code with `importlib.resources.files("plastid") / "data" / "<name>"`.
