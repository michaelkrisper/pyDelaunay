# pyDelaunay

Pure-Python Delaunay triangulation of scattered 3D points with linear Z interpolation.

[![CI](https://github.com/michaelkrisper/pyDelaunay/actions/workflows/ci.yml/badge.svg)](https://github.com/michaelkrisper/pyDelaunay/actions/workflows/ci.yml)
[![Version](https://img.shields.io/github/v/tag/michaelkrisper/pyDelaunay?sort=semver&label=version)](https://github.com/michaelkrisper/pyDelaunay/releases)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

`DelaunayMap` triangulates the (X, Y) projection of a point cloud using the
Bowyer-Watson algorithm and returns the Z value of any (X, Y) coordinate inside
the convex hull by interpolating on the plane of the enclosing triangle.

## Features

- Bowyer-Watson incremental triangulation, no dependencies outside the standard library
- Planar Z interpolation via `d[x, y]` (returns `None` outside the convex hull)
- Triangle vertex indices via `d.triangle_indices`, e.g. for plotting with matplotlib
- Single module (`delaunay.py`), easy to vendor

## Installation

```bash
pip install git+https://github.com/michaelkrisper/pyDelaunay.git
```

Or copy `delaunay.py` into your project.

## Usage

```python
from delaunay import DelaunayMap

d = DelaunayMap((0, 0, 0), (1, 0, 1), (0, 1, 1))

print(d[0.25, 0.25])            # 0.5
print(d[2, 2])                  # None (outside the convex hull)
print(list(d.triangle_indices)) # [(0, 1, 2)]
```

Points are `(x, y, z)` tuples; their (X, Y) coordinates must be unique and at
least three points are required.

### Demo

Running the module triangulates 200 random points (or the points of a CSV file
with `x,y,z` lines), prints timings and plots the mesh with matplotlib:

```bash
pip install matplotlib
python delaunay.py            # random points
python delaunay.py points.csv # your own points
```

## Development

```bash
pip install pytest ruff
pytest
ruff check .
```

CI runs both on every pull request. Every merge to `main` is tagged with the
next patch version; bundled releases with release notes come from the Release
workflow.

## Contributing

Changes start with an issue and go through a pull request; see
[CONTRIBUTING](CONTRIBUTING.md).

## License

MIT - see [LICENSE](LICENSE).
