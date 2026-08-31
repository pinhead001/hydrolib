# Contributing to HydroLib

HydroLib computes flood frequency estimates that people use to size structures and set
insurance rates. A wrong number here is worse than a missing one, so most of what follows is
about making errors visible rather than about style.

## Setup

```bash
git clone https://github.com/pinhead001/hydrolib
cd hydrolib
pip install -e ".[dev]"
```

That installs the **pinned** toolchain — black 24.10, isort 5.13, pytest 7.4, mypy 1.x. Install
the pins, not the latest: CI installs this same extra, and a bare `pip install black` resolves
black 26, which formats differently and will fail the build. This is not hypothetical; it went
red in August 2026 exactly this way.

Optionally, run the checks on every commit:

```bash
pip install pre-commit && pre-commit install
```

## The one command

```bash
make check       # lint + typecheck + test, in CI's order
```

Everything CI runs is a `make` target, so local and CI cannot drift:

| target | what it is |
|---|---|
| `make lint` | black and isort, check-only, over `hydrolib/ tests/ app/` |
| `make fmt` | the same two, applying changes |
| `make typecheck` | mypy over `hydrolib/` |
| `make test` | the suite, exactly as CI selects it |
| `make coverage` | the suite with coverage, enforcing the floor |
| `make smoke` | imports and exercises the Streamlit app (needs `-r app/requirements.txt`) |
| `make parity` | builds the Fortran extension and checks the golden files (needs gfortran + meson) |
| `make golden` | regenerates the parity golden files |

Two traps that have each cost a red build:

- `python -m pytest` puts the working directory on `sys.path`; the `pytest` console script that
  CI runs does not. Use `PYTHONSAFEPATH=1 python -m pytest`, which is what the Makefile does.
- CI tests Python 3.9–3.12. Green on one says nothing about the others — fixture plumbing in
  particular diverges (a `@staticmethod` fixture works on 3.11 and breaks collection on 3.9).

## Conventions

- Type hints on every signature; NumPy-format docstrings on public API.
- No bare `except:`. No `print()` in library code — use `logging.getLogger(__name__)`.
- Every public function gets tests: happy path plus at least one error case.
- All data in log10 space when interfacing with Fortran conventions.
- **Never edit anything under `vendor/`.** It is a verbatim copy of USGS peakfq 8.1.0, and it is
  the thing every parity comparison is made against. A change there does not fail a test; it
  silently redefines what "correct" means.

## Two ratchets

Both exist so a gate could be turned on today rather than left off indefinitely. Both are meant
to tighten. Neither is a place to park new debt.

**mypy** excludes six modules in `pyproject.toml`'s `[[tool.mypy.overrides]]` — 155 errors
between them. The other fourteen are gated. Fix a module, delete its line in the same PR. Do not
add a line: a module that is off the list stays off it.

**Coverage** has a floor of 66% (`COV_MIN` in the Makefile) against a measured 68%. Raise the
floor when you raise the coverage. `cli.py`, `plots.py` and `validation/reports.py` are still at
0%, and `hydrolib.batch` was too, which is how `analyze_sites` shipped raising `AttributeError`
for every site — a bare `except` turned it into `{"error": ...}` and nothing ever executed the
function to notice.

## Known-failing tests

A test for a defect you understand but have not fixed is `xfail(strict=True)`. Never skip one,
and never widen a tolerance to hide it. Strict means the build fails the moment it starts
passing, which is the alarm you want when someone accidentally fixes it. Say so in the PR when
you add one.

There are currently seven, all real and all documented in `TODO.md`: the weighted-skew and
confidence-interval defects, which both bottom out in the unported Fortran `var_mom`.

## Numerical changes

The EMA fixed point is ill-conditioned — a 1-ulp input perturbation moves converged outputs by
up to 3.2e-3, a condition number around 1e13. Two consequences:

- Do not "robustify" by rounding. Rounding perception thresholds to 12 decimals moved Big
  Sandy's at-site skew MSE by 2.2e-4.
- Do not verify a ported Fortran routine through `emafitpr` end to end. Through a fixed point
  that badly conditioned, a correct routine can look wrong and a wrong one can look right.
  `tests/fortran_parity/test_fortran_oracles.py` calls `mseg_all_sub`, `detratsub`, `var_mom`
  and `moms_p3` directly; add your routine there.

If you change anything under `build_fortran/` or regenerate goldens, run `make parity` and say
in the PR what moved and by how much.

## Pull requests

Say what changed and how you know it is right. For a numerical change that means before/after
numbers against the reference, not "tests pass". CI must be green: lint, mypy, the 3.9–3.12
matrix, coverage, the Streamlit smoke test, and Fortran parity.

## Releases

Maintainers only, and tag-driven:

```bash
bump2version patch      # or minor / major
git push && git push --tags
```

`.github/workflows/release.yml` refuses to release unless the tag, `pyproject.toml` and
`.bumpversion.cfg` all agree on the version, then builds, tests, checks the wheel actually
contains `py.typed` and the packaged benchmark data, and publishes a GitHub release. PyPI
publishing is opt-in and stays off until a `PUBLISH_TO_PYPI` repository variable is set to
`true` and trusted publishing is configured.
