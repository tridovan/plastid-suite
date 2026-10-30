# Cutting a release (M9.6)

**Why this file exists.** A release is the one operation in this project that cannot be undone: a
DOI is minted per version, and a PyPI filename is taken the moment it is uploaded. So the steps are
written down in order, with the check that has to pass between them, and with the two things that
are **decisions, not steps** — the version number and who is named as an author.

## Before anything is published

Publication is deferred while the [known correctness issues](known_issues.md) remain unresolved
for the advertised scope. Passing the existing suite alone is not product acceptance: each fixed
defect needs an independently justified regression and a recorded before/after result. Record
missing tools and panels; skipped acceptance cases remain unverified.
Use the [pipeline coverage and acceptance criteria](pipeline.md) to check reference preparation,
both genome entry paths, analyses, figures/exports/QC and submission; do not infer end-to-end
coverage merely because each component has some tests. Deferred branches are not release promises.

| Check | Command | What it proves |
|---|---|---|
| The suite | `make check`, with required tools present and skips reviewed | lint, fast tests and golden cases that actually ran; this does not replace input/output acceptance |
| The panel | `make acceptance`, on a machine that holds `.panel/` | the tests marked `panel` ran on the validation panel with no skip — the tier CI cannot run |
| The package | `uv build` then install the wheel into a **fresh** venv | the artifact works somewhere that has nothing of the dev's in it |
| The container | `docker build -f docker/Dockerfile -t plastid-suite .` | the golden numbers reproduce **inside** the image — one input, one set of numbers, three environments (WSL · CI · Docker) |
| The citation | `uv tool run cffconvert --validate` | `CITATION.cff` parses against CFF 1.2.0 |

## The steps, in order

1. **Decide the version**. It goes in `pyproject.toml`, `README.md`'s status line and
   `CITATION.cff`. Check all three agree — `tests/unit/test_product_docs.py` fails if the README
   status line does not name the `pyproject.toml` version, and its `test_the_version_check_reports_a_mismatch`
   runs the same predicate against a stale and an absent line, so the check is known to fail when
   it should. The development-status classifier is a separate decision from the number: it is
   `4 - Beta` while the correctness issues are open, and the README says so.
2. **Decide the authors** — `CITATION.cff` currently names "PlastidToolSuite contributors"
   deliberately. A Zenodo deposit records its authors permanently, so the real list is fixed
   **before** step 5, by the team.
3. **Commit the release state** on `main`, then
   `git tag v<version> && git push origin main --tags`. `main` is one commit per phase — a release
   tag points at a state somebody can check out and trust.
4. **Build the artifacts**: `uv build` → `dist/*.whl`, `dist/*.tar.gz`. Keep them: the sdist is
   what a reviewer reads, and both are what the checks above were run against.
5. **Publish to PyPI** (the dev, explicitly — outward-facing and effectively permanent):

   ```bash
   uv publish            # or: twine upload dist/*
   ```

   The name `plastid-suite` on PyPI is claimed by this upload; check it is still free first.
6. **GitHub release** for the tag, with the artifacts attached. Zenodo's GitHub integration
   archives it and mints the DOI. **Note:** Zenodo's integration works on public repositories —
   this repository is private, so either the release is cut after the repository is made public, or
   the tarball is uploaded to Zenodo by hand. *Not yet done; decided by the dev at release time.*
7. **Put the DOI back into the tree** — `CITATION.cff` (`doi:`) and `README.md` —
   and commit. A DOI that lives only in a web page is a DOI the next reader has to go looking for.
8. **The how-to video**: recorded by the team, linked from the README. No code produces it, and the
   release is not "done" without it (M9.6).

## What must NOT happen

* **Publishing a prebuilt Docker image.** It contains MISA, whose licence does not permit
  redistribution (`docs/licensing.md`, D-q). The project ships the *recipe*; each user builds their
  own copy.
* **A version bump without the README and CITATION.cff**, or a tag without the artifacts built from
  the same commit — the three must describe one state.
* **Publishing from a working tree that is not clean**, or from a branch other than `main`.
