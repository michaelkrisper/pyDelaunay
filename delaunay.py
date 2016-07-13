#! /usr/bin/env python3
# coding=utf-8
"""
An implementation of delaunay triangulation.
"""
from collections import Counter, namedtuple
from math import isclose
import time
from typing import Tuple, Iterable, List

__author__ = "Michael Krisper"
__email__ = "michael.krisper@gmail.com"
__date__ = "2016-07-12"
__python_version__ = "3.5"

Point = namedtuple("Point", "X Y Z index")


class Point2:
    def __init__(self, x: float, y: float, z: float, index: int):
        self.X = float(x)
        self.Y = float(y)
        self.Z = float(z)
        self.index = index

    def __repr__(self):
        return "({X:.7}, {Y:.7}, {Z:.7})".format_map(vars(self))

    def __eq__(self, other):
        return isclose(self.X, other.X) and isclose(self.Y, other.Y) and isclose(self.Z, other.Z)

    def __hash__(self):
        return hash((self.X, self.Y, self.Z))

    def __lt__(self, other):
        return (self.X, self.Y, self.Z).__lt__((other.X, other.Y, other.Z))

class Edge:
    def __init__(self, p1: Point, p2: Point):
        self.P1 = p1
        self.P2 = p2

    def __hash__(self):
        return hash((min(self.P1, self.P2), max(self.P1, self.P2)))

    def __eq__(self, other):
        return (self.P1, self.P2) == (other.P1, other.P2) or (self.P1, self.P2) == (other.P2, other.P1)

    def __repr__(self):
        return "({P1} -> {P2})".format_map(vars(self))


class Triangle:
    def __init__(self, p1: Point, p2: Point, p3: Point):
        if p1.X == p2.X == p3.X or p1.Y == p2.Y == p3.Y:
            raise Exception(
                "triangle is coplanar --> cannot extrapolated by a plane.\np1:{p1}\np2:{p2}\np3:{p3}".format_map(
                    vars()))
        self.P1 = p1
        self.P2 = p2
        self.P3 = p3

        ab = Point(p2.X - p1.X, p2.Y - p1.Y, p2.Z - p1.Z, -1)
        ac = Point(p3.X - p1.X, p3.Y - p1.Y, p3.Z - p1.Z, -1)
        self.X = ab.Y * ac.Z - ab.Z * ac.Y
        self.Y = ab.Z * ac.X - ab.X * ac.Z
        self.Z = ab.X * ac.Y - ab.Y * ac.X
        self.W = p1.X * self.X + p1.Y * self.Y + p1.Z * self.Z

    def __repr__(self):
        return "<{P1}, {P2}, {P3}>".format_map(vars(self))

    def shares_vertex(self, tr: 'Triangle') -> bool:
        return {self.P1, self.P2, self.P3} & {tr.P1, tr.P2, tr.P3}

    def contains_in_circumcircle(self, p: Point) -> bool:
        p1X, p1Y = self.P1.X, self.P1.Y
        p2X, p2Y = self.P2.X, self.P2.Y
        p3X, p3Y = self.P3.X, self.P3.Y
        p4X, p4Y = p.X, p.Y
        p4X2, p4Y2 = p4X ** 2, p4Y ** 2
        a, b, c = p1X - p4X, p1Y - p4Y, (p1X ** 2 - p4X2) + (p1Y ** 2 - p4Y2)
        d, e, f = p2X - p4X, p2Y - p4Y, (p2X ** 2 - p4X2) + (p2Y ** 2 - p4Y2)
        g, h, i = p3X - p4X, p3Y - p4Y, (p3X ** 2 - p4X2) + (p3Y ** 2 - p4Y2)
        det = a * (e * i - f * h) + b * (f * g - d * i) + c * (d * h - e * g)
        return det > 0 and not isclose(d, 0)

    def is_inside(self, x: float, y: float) -> bool:
        v0X, v0Y = self.P3.X - self.P1.X, self.P3.Y - self.P1.Y
        v1X, v1Y = self.P2.X - self.P1.X, self.P2.Y - self.P1.Y
        v2X, v2Y = x - self.P1.X, y - self.P1.Y
        dot00 = v0X ** 2 + v0Y ** 2
        dot01 = v0X * v1X + v0Y * v1Y
        dot02 = v0X * v2X + v0Y * v2Y
        dot11 = v1X ** 2 + v1Y ** 2
        dot12 = v1X * v2X + v1Y * v2Y
        inv_denom = 1.0 / (dot00 * dot11 - dot01 ** 2)
        u = (dot11 * dot02 - dot01 * dot12) * inv_denom
        v = (dot00 * dot12 - dot01 * dot02) * inv_denom
        return (u > 0 or isclose(u, 0)) and (v > 0 or isclose(v, 0)) and (u + v < 1 or isclose(u + v, 1))

    def interpolate(self, x: float, y: float) -> float:
        return (self.W - self.X * x - self.Y * y) / self.Z


class DelaunayMap:
    def __init__(self, *_points: Iterable[Tuple[float, float, float]]):
        points = list(Point(p[0], p[1], p[2], i) for i, p in enumerate(_points))
        if len(points) != len({(p.X, p.Y) for p in points}):
            raise Exception("(X,Y) coordinates of the points must be unique!")
        if len(points) < 3:
            raise Exception("There must be at least 3 points!")

        min_x, min_y = min(p.X for p in points), min(p.Y for p in points)
        max_x, max_y = max(p.X for p in points), max(p.Y for p in points)
        supertriangle = Triangle(Point(min_x - 1, min_y - 1, 0, -1), Point(max_x * 3, min_y - 1, 0, -1),
                                 Point(min_x - 1, max_y * 3, 0, -1))
        triangles = [supertriangle]  # type: List[Triangle]
        point_count = 0
        for p in points:
            point_count += 1
            container_triangles = [tr for tr in triangles if tr.contains_in_circumcircle(p)]  # type: List[Edge]
            edges = []  # type: List[Edge]
            for tr in container_triangles:  # type: Triangle
                triangles.remove(tr)
                edges += [Edge(tr.P1, tr.P2), Edge(tr.P2, tr.P3), Edge(tr.P3, tr.P1)]
            unique_edges = (edge for edge, count in Counter(edges).items() if count == 1)
            new_triangles = (Triangle(e.P1, e.P2, p) for e in unique_edges)
            triangles += new_triangles

            if len(triangles) != 2 * (point_count + 3) - 2 - 3:
                raise Exception(
                    "Triangulation invariant violated! Triangle count and point count doesn't fit together.")

        self._triangles = [tr for tr in triangles if not tr.shares_vertex(supertriangle)]  # type: List[Triangle]

    def __getitem__(self, xy: Tuple[float, float]) -> float:
        x, y = xy
        for tr in self._triangles:
            if tr.is_inside(x, y):
                return tr.interpolate(x, y)

    @property
    def triangles_indizes(self):
        for tr in self._triangles:  # type: Triangle
            yield tr.P1.index, tr.P2.index, tr.P3.index


if __name__ == '__main__':
    with open("map.csv") as f:
        points = [tuple(float(x) for x in line.strip().split(",")) for line in f.readlines()]

    # random.shuffle(points)
    start = time.clock()
    d = DelaunayMap(*points)
    end = time.clock()
    print(end - start)

    start = time.clock()
    for x, y, z in points:
        assert z == d[x, y]
    end = time.clock()
    print(end - start)

    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    import numpy as np

    arr = np.array(points).astype("d")
    x, y, _ = arr.transpose()
    triang = mtri.Triangulation(x, y, list(d.triangles_indizes))
    xi, yi = np.meshgrid(np.linspace(0, 5, 200), np.linspace(0, 5, 200))

    plt.triplot(triang, 'ko-')

    plt.show()
