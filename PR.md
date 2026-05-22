## Title
Modernize packaging/CI and overhaul docs for RTD rendering and CIGAR guidance

## Summary
- Modernizes project packaging and metadata by moving to `pyproject.toml`, removing legacy `setup.py`, and updating development requirements.
- Adds/updates CI and documentation infrastructure (`.github/workflows/ci.yml`, `.readthedocs.yaml`, Sphinx config/requirements) for current Python and Read the Docs workflows.
- Expands user-facing documentation substantially, including a practical `usage` guide and a new dedicated CIGAR format explainer page.
- Normalizes API docstring examples across core modules so fixed-width alignment examples render correctly and consistently on RTD.
- Updates tox matrix to current supported Python versions (`3.8`-`3.12`) to match project classifiers.
- Applies a full repository flake8 cleanup across package and tests, addressing style debt and unused-symbol issues.

## Commits Included (main...publish-prep)
1. `864c0f2` modernizing
2. `4418b42` fixing docs for fixed width setup
3. `ed10800` fixes for sphinx warnings
4. `60fe528` unpinning
5. `a6b88ea` pointing at correct repo
6. `dba9d6f` adding cigar explanation and basic usage
7. `023158b` updating old tox
8. `858c865` flake8 massive fix

## Notable Changes by Area

### Packaging and metadata
- Added `pyproject.toml` with project metadata and dependency groups.
- Removed legacy `setup.py`.
- Consolidated dev tooling into the `dev` extra in `pyproject.toml` (retired `requirements_dev.txt`).
- Updated package version metadata reference points in `cigarmath/__init__.py`.

### CI and release/doc infrastructure
- Added GitHub Actions workflow at `.github/workflows/ci.yml`.
- Added Read the Docs config at `.readthedocs.yaml`.
- Updated `tox.ini` env matrix to `py38, py39, py310, py311, py312, flake8`.

### Documentation structure and quality
- Expanded docs pages: `docs/api.rst`, `docs/installation.rst`, `docs/readme.rst`, `docs/history.rst`, `docs/index.rst`.
- Added new CIGAR reference page: `docs/cigar.rst`.
- Rewrote `docs/usage.rst` into a workflow-based guide with runnable examples and expected outputs.
- Updated docs build requirements: `docs/requirements.txt`.
- Tuned Sphinx config in `docs/conf.py` for stable rendering and warning cleanup.

### API docstring rendering fixes
- Updated fixed-width examples to explicit literal-block style in:
  - `cigarmath/block.py`
  - `cigarmath/cigarmath.py`
  - `cigarmath/clipping.py`
  - `cigarmath/inference.py`
  - `cigarmath/mapping.py`
  - `cigarmath/iterators.py`
  - `cigarmath/combine.py`
  - `cigarmath/pileup.py`
  - `cigarmath/conversions.py`

### Lint and style sweep
- Ran `flake8` cleanup across `cigarmath/` and `tests/`.
- Resolved unused-import/redefinition issues in public exports (`cigarmath/__init__.py`) and module imports.
- Removed checkpoint artifact from lint surface (`.ipynb_checkpoints` excluded / cleaned).
- Standardized formatting and whitespace across tests and core modules.

## Diff Scope
- **42 files changed**
- **1476 insertions**, **831 deletions**

## Test Plan
- [x] Build docs locally:
  - `python -m sphinx -b html docs docs/_build/html`
- [x] Validate no Sphinx build warnings in latest pass.
- [x] Run style checks:
  - `flake8 cigarmath tests`
- [ ] Run tox across available local interpreters:
  - `tox -e py38,py39,py310,py311,py312,flake8`
- [ ] Run unit tests directly:
  - `pytest`

## Risk / Reviewer Notes
- Largest review surface is docs, packaging modernization, and the broad lint-driven formatting sweep.
- Core behavioral intent should be unchanged; most additional churn is style/format and import hygiene.
- Recommend reviewers focus on:
  - Packaging metadata correctness (`pyproject.toml`)
  - CI/tox Python version alignment
  - RTD/Sphinx rendering behavior for fixed-width examples
  - Spot-checking representative tests/modules touched by flake8 normalization
