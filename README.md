# pyDelaunay

[![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A basic Delaunay triangulation implementation in Python.

This project provides a simple implementation of the Delaunay triangulation algorithm for a set of 2D points (with Z values). It can be used to interpolate the Z value for any point within the convex hull of the input points.

## Usage

The main class is `DelaunayMap`, which is initialized with a list of 3D points. You can then get the interpolated Z value for any (X, Y) coordinate by using the `__getitem__` method.

Here is an example:

```python
from delaunay import DelaunayMap

# Create a DelaunayMap with three points
d = DelaunayMap((0, 0, 0), (1, 0, 1), (0, 1, 1))

# Interpolate the Z value at (0.25, 0.25)
z = d[0.25, 0.25]

print(z)  # Output: 0.5
```

## Running the Tests

To run the tests, simply execute the `test.py` file:

```bash
python test.py
```

This will run a series of unit tests to ensure that the triangulation and interpolation are working correctly.

## Demonstration

The `delaunay.py` file also includes a demonstration that reads points from a `map.csv` file and plots the triangulation using `matplotlib`. To run the demonstration, you will need to have `matplotlib` installed and a `map.csv` file in the same directory.

```bash
pip install matplotlib
python delaunay.py
```

## Requirements

- Python 3.9 or newer
- `matplotlib` — only for the optional plotting demonstration

## Project structure

| File | Purpose |
|---|---|
| `delaunay.py` | `DelaunayMap` — triangulation and Z interpolation |
| `test.py` | Unit tests for triangulation and interpolation |
| `setup.py` | Packaging metadata |

## License

MIT — see [LICENSE](LICENSE).
