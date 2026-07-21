# GeoEq Release Guide

Step-by-step checklist for publishing a new GeoEq release to PyPI,
GitHub, and the documentation sites.

GeoEq follows [Semantic Versioning](https://semver.org/):

| Bump | When |
| :-- | :-- |
| **Patch** (`0.1.3` → `0.1.4`) | Bug fixes, docs, metadata — no API change |
| **Minor** (`0.1.x` → `0.2.0`) | New functions or modules, backwards-compatible |
| **Major** (`0.x` → `1.0.0`) | Breaking API changes |

---

## 1. Pre-flight checks

```bash
cd geoeq

# Working tree must be clean and synced
git status
git pull origin main

# Full test suite must pass
pytest                        # → all tests green
```

Do not proceed if any test fails.

## 2. Bump the version

The version lives in **three** places — update all of them:

| File | Field |
| :-- | :-- |
| `pyproject.toml` | `version = "X.Y.Z"` |
| `src/geoeq/__init__.py` | `__version__ = "X.Y.Z"` |
| `README.md` | Capabilities intro + citation BibTeX `version = {X.Y.Z}` |

Quick check that nothing was missed:

```bash
grep -rn "OLD_VERSION" --include="*.py" --include="*.toml" --include="*.md" . \
  | grep -v CHANGELOG | grep -v .git
```

## 3. Update the CHANGELOG

Add a new section at the top of `CHANGELOG.md` following
[Keep a Changelog](https://keepachangelog.com/):

```markdown
## [X.Y.Z] - YYYY-MM-DD

### Added / Changed / Fixed / Removed
- One bullet per user-visible change.

### Note
State test count and whether any API changed.
```

## 4. Build the distribution

```bash
# Clean old artefacts, then build sdist + wheel
rm -rf dist/ build/
python -m build

# Validate metadata renders correctly on PyPI
python -m twine check dist/*      # → PASSED, PASSED
```

Verify the licence files ship inside the sdist:

```bash
tar -tzf dist/geoeq-X.Y.Z.tar.gz | grep -E "LICENSE|NOTICE"
```

## 5. (Optional) Rehearse on TestPyPI

For risky releases, upload to TestPyPI first:

```bash
python -m twine upload --repository testpypi dist/*
pip install --index-url https://test.pypi.org/simple/ \
            --extra-index-url https://pypi.org/simple/ geoeq==X.Y.Z
python -c "import geoeq as ge; print(ge.__version__)"
```

## 6. Upload to PyPI

```bash
python -m twine upload dist/*
```

Username is `__token__`; the password is your PyPI API token
(scoped to the `geoeq` project). Never commit the token anywhere.

## 7. Commit, tag, and push

```bash
git add pyproject.toml src/geoeq/__init__.py README.md CHANGELOG.md
git commit -m "chore(release): vX.Y.Z — <one-line summary>"
git push origin main

git tag vX.Y.Z
git push origin vX.Y.Z
```

## 8. Create the GitHub release

```bash
gh release create vX.Y.Z \
  --title "GeoEq vX.Y.Z" \
  --notes-file <(sed -n '/## \[X.Y.Z\]/,/^---$/p' CHANGELOG.md) \
  dist/geoeq-X.Y.Z.tar.gz dist/geoeq-X.Y.Z-py3-none-any.whl
```

Or draft it in the GitHub UI — paste the CHANGELOG section as the
release notes and attach both `dist/` files.

## 9. Update the websites

Version numbers and badges appear on more than the package:

| Site | Repo | What to update |
| :-- | :-- | :-- |
| geoeq.org (landing) | `geoeqlanding` | Hero version badge, meta tags |
| docs.geoeq.org | `geoeqdocs` | `installation.md` version pins |
| geoeq.org (Pages) | `geoeq.github.io` | Stat cards, user guide |

The PyPI badges in the README refresh automatically.

## 10. Post-release verification

Wait a few minutes for PyPI to propagate, then verify from a clean
environment:

```bash
python -m venv /tmp/geoeq-verify && source /tmp/geoeq-verify/bin/activate
pip install geoeq
python -c "
import geoeq as ge
print(ge.__version__)                          # → X.Y.Z
print(ge.bearing_capacity(c=10, gamma=18, Df=1, B=2, phi=30)['q_u'])
"
deactivate && rm -rf /tmp/geoeq-verify
```

Finally check the [PyPI project page](https://pypi.org/project/geoeq/)
renders the README, licence, and links correctly.

---

## Quick reference — one release, ten commands

```bash
pytest                                   # 1  all green
# bump version in 3 files + CHANGELOG    # 2-3
rm -rf dist/ build/ && python -m build   # 4
python -m twine check dist/*             # 4
python -m twine upload dist/*            # 6
git add -A && git commit -m "chore(release): vX.Y.Z"
git push origin main                     # 7
git tag vX.Y.Z && git push origin vX.Y.Z # 7
gh release create vX.Y.Z --title "GeoEq vX.Y.Z"  # 8
pip install geoeq && python -c "import geoeq; print(geoeq.__version__)"  # 10
```

---

*GeoEq — geotechnical engineering, solved in Python.
[geoeq.org](https://geoeq.org) · [PyPI](https://pypi.org/project/geoeq/) ·
Apache 2.0*
