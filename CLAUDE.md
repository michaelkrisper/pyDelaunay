# CLAUDE.md

Pure-Python Delaunay triangulation (Bowyer-Watson) of the (X, Y) projection of 3D
points with planar Z interpolation. One module, `delaunay.py`, Python 3.9+, standard
library only; matplotlib/numpy only for the `__main__` demo. Installable via
`pyproject.toml` (setuptools, `py-modules = ["delaunay"]`); nothing is deployed.

## Commands

```bash
pip install pytest ruff
pytest                           # test_delaunay.py, < 1 s
ruff check .                     # config in pyproject.toml (E, F, W, I, B; 120 cols)
pip install matplotlib && python delaunay.py [points.csv]   # demo plot, x,y,z per line
```

CI (`.github/workflows/ci.yml`) runs ruff on 3.14 and pytest on 3.9 and 3.14.

## Architecture

- `DelaunayMap(*points)`: validates input (>= 3 points, unique (X, Y)), builds a
  supertriangle, inserts points one by one (circumcircle test, re-triangulate the hole
  from edges seen exactly once), then drops triangles touching the supertriangle.
- `d[x, y]` returns the interpolated Z of the first triangle whose projection contains
  the point, `None` outside the convex hull. `triangle_indices` yields input indices.
- `Point` (namedtuple, carries the input `index`) is what the algorithm uses;
  supertriangle vertices have index `-1`. `Point2` is unused legacy.

## Rules and pitfalls

- Keep the supertriangle at `size * 100` around the bounding box. A smaller one lost
  convex hull triangles (a three-point map returned `None`).
- `contains_in_circumcircle()` multiplies by the triangle orientation, so it works for
  clockwise and counter-clockwise triangles; don't assume a winding order.
- The triangle-count invariant check after each insertion is deliberate; if it fires,
  fix the algorithm instead of removing the check.
- `triangles_indizes` is a backwards-compatible alias of `triangle_indices`; keep it.
- Tests check properties on random point sets (Euler triangle count, empty
  circumcircles, every point used, a plane reproduced exactly); add a regression test
  for every geometry fix.
- Stay on the 3.9 syntax level (`typing.List` etc. are fine, no `X | Y` unions) and
  standard library only.
- `version` in `pyproject.toml` is static; releases are git tags only.

## Workflow

- Every change starts with a GitHub issue. Branch `<issue>-<slug>` from `main`,
  open a PR with `Closes #<issue>`. Never push to `main` directly.
- PRs auto-merge (squash) once CI and the "Linked ticket" check pass; drafts are not merged.
- Versions exist only as git tags `vX.Y.Z`. Every merge to `main` tags the next patch;
  never edit version numbers in files. Minor/major releases: Actions -> Release.
- Nothing is deployed; a merge to `main` only creates the patch tag.
- Details: CONTRIBUTING.md.
