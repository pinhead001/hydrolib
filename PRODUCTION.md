# HydroLib Production Readiness Status

**Status**: ✅ **Phase 1 and the Phase 2 tooling complete** | **Next Phase**: docs, license headers, integration tests

Last updated: 2026-08-31

---

## Phase 1: Minimum Viable Changes (MVC) — ✅ COMPLETE

All critical governance and security measures are in place.

### Completed

| Item | Deliverable | Status |
|------|-------------|--------|
| **Remove setup.py** | Deleted; pyproject.toml is single source of truth | ✅ |
| **Repair `bump2version`** | Config said 0.0.3 while every file said 0.2.0, so it raised `VersionNotFoundException` and no release could be cut | ✅ |
| **PEP 561 marker** | `hydrolib/py.typed` added and packaged | ✅ |
| **Code ownership** | `.github/CODEOWNERS` — all code requires review | ✅ |
| **Security policy** | `SECURITY.md` — vulnerability reporting workflow | ✅ |
| **Branch protection** | `main` requires: 1 PR review + all CI checks pass | ✅ |
| **Author metadata** | `pyproject.toml` updated with author, classifiers, keywords | ✅ |
| **Dependency pinning** | Dev tools locked: black <25, pytest <8, isort <6, flake8 <7, pytest-cov <5 | ✅ |
| **Changelog** | `CHANGELOG.md` (Keep a Changelog format) with v0.1.0–v0.2.0 history | ✅ |

### Commit
- **Hash**: `77cc03b9556229a38e5e2d30770b5f8b93ef4eb5`
- **Message**: `chore: production readiness MVC - delete setup.py, add governance files, pin dependencies`
- **Impact**: 4 files pushed to main

---

## Phase 2: Medium-Priority Enhancements (Recommended: 2–3 weeks)

These should be completed before widespread production adoption.

### Type Checking & Linting — ✅ mypy and pre-commit done

| Task | Priority | Effort | Status |
|------|----------|--------|--------|
| ~~Add `mypy` to CI and enforce type hints~~ | HIGH | est. 2–3 days; took ~1 hour | ✅ `typecheck` job, `make typecheck` |
| ~~Create `.pre-commit-config.yaml` for local enforcement~~ | MEDIUM | est. 1 day | ✅ pinned to the dev extra's versions |
| ~~Add type hints to all function signatures~~ — already done (99%) | — | — | ✅ |
| Add `pydocstyle` to CI for docstring coverage | MEDIUM | 1–2 days | open |

**What mypy found on its first run**, which is the answer to "highest ROI": 172 errors, of
which four were real defects in shipped code — `historical_peaks`/`perception_thresholds` on
`Bulletin17C` and both date arguments to `USGSgage.download_daily_flow()` declared
non-optional with a `None` default, so a type checker rejected the documented call; and two
`-> "pd.DataFrame"` annotations naming a module never imported, unresolvable by mypy, an IDE
or `typing.get_type_hints()`.

**Fourteen of twenty modules are gated.** The other six — `__init__`, `bulletin17c`,
`freq_plot`, `hydrograph`, `lowflow`, `report` — hold 155 errors and are excluded in
`pyproject.toml`'s `[[tool.mypy.overrides]]`. That list is a ratchet: fix a module, delete
its line. Nothing may be added to it.

---

### Release Process & PyPI Publishing — ✅ workflow done, publishing opt-in

| Task | Priority | Effort | Status |
|------|----------|--------|--------|
| ~~Create `.github/workflows/release.yml`~~ | HIGH | 1 day | ✅ tag-driven |
| ~~Document version bump workflow (using `bump2version`)~~ | MEDIUM | 2 hours | ✅ `CONTRIBUTING.md` |
| ~~Add `PYPI_TOKEN` secret~~ — superseded | HIGH | — | ✅ n/a: trusted publishing, no long-lived token |
| Enable PyPI publishing | HIGH | 30 min | open — needs the owner |
| Test the workflow against TestPyPI | MEDIUM | 1 hour | open |

The workflow builds, tests, and publishes a GitHub release on a `v*` tag. Two gates run
before anything is published: the tag, `pyproject.toml` and `.bumpversion.cfg` must agree on
the version — the check that would have caught the 0.0.3-vs-0.2.0 drift above — and the built
wheel must actually contain `py.typed` and the packaged benchmark data, both of which have
gone missing before.

PyPI publishing is deliberately **off** until someone sets the `PUBLISH_TO_PYPI` repository
variable to `true` and registers hydrolib as a trusted publisher, so a tag push cannot
publish by accident. Trusted publishing is why there is no `PYPI_TOKEN`: OIDC issues a
short-lived credential per run, leaving no long-lived secret to leak.

**Not verified end to end.** Tag pushes return HTTP 403 from the agent environment, so the
workflow has never run. What *was* verified locally, with the pinned toolchain: the version
guard accepts `v0.2.0` and rejects `v0.3.0`; `python -m build` produces both artifacts;
`twine check` passes on each; and the wheel-contents assertion passes on the real wheel.

---

### Test Coverage Reporting — ✅ done, floor at 66%

| Task | Priority | Effort | Status |
|------|----------|--------|--------|
| ~~Add `pytest-cov` to a CI job~~ | HIGH | 30 min | ✅ `coverage` job, HTML uploaded as an artifact |
| ~~Set a coverage threshold in CI~~ | MEDIUM | 30 min | ✅ 66% floor, `make coverage` |
| Add a badge to README.md | LOW | 15 min | open |
| Integrate with codecov.io (optional) | LOW | 1 hour | open |

**The threshold is 66%, not the 80% this document proposed.** Measured coverage is 67.87%;
80% would have failed CI on the first run. A floor has to be under the floor. Raise it as
coverage rises — `COV_MIN` in the Makefile — and verify it in both directions: at 66 the
gate passes, at 70 it fails with a non-zero exit.

Coverage went 65% → 68% in the course of adding the gate, because measuring it turned up
`hydrolib.batch` at 0% and `analyze_sites()` broken for every input (see Fixed in the
changelog). Still at 0%: `cli.py`, `plots.py`, `validation/reports.py`.

---

### Documentation & License Headers

| Task | Priority | Effort | Impact |
|------|----------|--------|--------|
| Add MIT license headers to all `.py` files | MEDIUM | 2 hours | Legal clarity; industry standard |
| ~~Add `CONTRIBUTING.md` with development workflow~~ | MEDIUM | 2 hours | ✅ at the repo root, where GitHub links it |
| Expand API documentation with examples per module | MEDIUM | 4–6 hours | Better discoverability |
| Add deployment guide for Streamlit Cloud | LOW | 2 hours | Users can self-host easily |

---

### Pre-Commit Hooks — ✅ config added

| Task | Priority | Effort | Status |
|------|----------|--------|--------|
| ~~Create `.pre-commit-config.yaml`~~ | LOW | 1 day | ✅ black/isort/mypy + file hygiene |
| ~~Document `pre-commit install`~~ | LOW | 30 min | ✅ `CONTRIBUTING.md` |

The black and isort revs must match the pins in the dev extra; the config says so, because
installing different versions is precisely how CI lint went red in August 2026. mypy runs as
a `local` hook against your own environment rather than through `mirrors-mypy`, which would
resolve an isolated one and disagree with CI. `vendor/` is excluded from the whitespace
hooks — it is a verbatim reference copy, and a stray newline there silently invalidates every
parity comparison.

**Not verified**: `pre-commit run --all-files` has not been executed here (the hook repos
need network fetches from the sandbox). The YAML parses and the local mypy hook runs the same
command as `make typecheck`, which is verified.

---

## Phase 3: Production Hardening (4–6 weeks)

Advanced safety measures for high-stakes deployments.

| Task | Priority | Effort | Impact |
|------|----------|--------|--------|
| **Automated security scanning** | MEDIUM | 2 days | Detect known vulnerabilities in dependencies |
| ~~Add `dependabot` config (`.github/dependabot.yml`)~~ | MEDIUM | 2 hours | ✅ pip + actions, monthly, linters grouped |
| Add SBOM (Software Bill of Materials) generation | LOW | 1 day | Supply chain transparency |
| Add integration tests (NWIS API calls in staging) | HIGH | 5–7 days | Catch breaking API changes early |
| Performance benchmarking CI job | LOW | 3–4 days | Detect regressions in computation time |
| Docker container for reproducible deployments | LOW | 2–3 days | Simplifies Streamlit Cloud deployment |

---

## Current Production Gaps & Risk Assessment

### 🟡 Medium Risk (Can go to production with caution)

1. **No integration tests with live USGS NWIS API**
   - Currently marked `requires_network` and skipped in CI
   - **Risk**: Breaking API changes from USGS will only surface in production
   - **Mitigation**: Add staging environment tests; subscribe to USGS API notifications
   - **Timeline**: 1–2 weeks

2. **Type hints present but unenforced**
   - Corrected 2026-08-31: the earlier claim that this codebase has no type
     hints was wrong. 210 of 213 functions in `hydrolib/` carry annotations
     (99%), as `CLAUDE.md` requires. The real gap was that `hydrolib/py.typed`
     was missing, so downstream type checkers ignored all of them — a
     one-file fix, now made.
   - What remains is *enforcement*: nothing checks the annotations are correct
     or that new code keeps them.
   - **Mitigation**: ✅ mypy is in CI. Much cheaper than the 3–5 days
     estimated below, because the annotations already exist.
   - **Timeline**: hours, not days

3. **No automated dependency updates**
   - Manual responsibility to update numpy, pandas, scipy, etc.
   - **Risk**: Missing security patches, using outdated versions
   - **Mitigation**: Enable Dependabot (Phase 3)
   - **Timeline**: 2 hours

4. **Limited Streamlit app testing**
   - `app/` is tested but with mocked/smoke tests only
   - No end-to-end UI testing
   - **Risk**: Streamlit updates may break UI
   - **Mitigation**: Document UI test checklist; add visual regression testing (optional)
   - **Timeline**: 2–3 days (optional)

### 🟢 Low Risk (Already addressed)

- ✅ **Code governance**: CODEOWNERS in place
- ✅ **Security reporting**: SECURITY.md defined
- ✅ **Dependency stability**: Dev tools pinned
- ✅ **Branch protection**: No direct commits to main; PRs require review + CI
- ✅ **Reproducibility**: Makefile ensures local = CI
- ✅ **Version history**: CHANGELOG.md maintained

---

## Recommended Rollout Timeline

### **Immediate (Next 24–48 hours)**
- ✅ Merge MVC changes (branch protection now active)
- Document rollout plan in team communications
- Brief team on new PR review requirements

### **Week 1–2 (Medium-priority)**
- ~~Add `mypy` type checking~~ ✅
- ~~Add test coverage enforcement~~ ✅
- ~~Set up the release workflow~~ ✅ (publishing to PyPI still needs enabling)
- Test the release process against TestPyPI

### **Week 3 (Nice-to-have)**
- Add Dependabot for dependency updates
- Add license headers
- Expand documentation

### **Before 1.0 Release (Milestone)**
- Add integration tests with live NWIS API (in staging)
- Achieve 80%+ test coverage
- Document performance SLAs
- Security audit (optional: third-party review)

---

## Files Changed (MVC Phase)
