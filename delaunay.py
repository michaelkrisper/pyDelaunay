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

Point = namedtuple("Point", "X Y Z")


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

    def __repr__(self):
        return "<{P1}, {P2}, {P3}>".format_map(vars(self))

    def contains_in_circumcircle(self, p: Point) -> bool:
        p1 = self.P1
        p2 = self.P2
        p3 = self.P3
        p4 = p
        a, b, c = p1.X - p4.X, p1.Y - p4.Y, (p1.X ** 2 - p4.X ** 2) + (p1.Y ** 2 - p4.Y ** 2)
        d, e, f = p2.X - p4.X, p2.Y - p4.Y, (p2.X ** 2 - p4.X ** 2) + (p2.Y ** 2 - p4.Y ** 2)
        g, h, i = p3.X - p4.X, p3.Y - p4.Y, (p3.X ** 2 - p4.X ** 2) + (p3.Y ** 2 - p4.Y ** 2)
        det = a * (e * i - f * h) + b * (f * g - d * i) + c * (d * h - e * g)
        return det > 0 and not isclose(d, 0)

    def shares_vertex(self, tr: 'Triangle') -> bool:
        return {self.P1, self.P2, self.P3} & {tr.P1, tr.P2, tr.P3}

    def is_inside(self, x: float, y: float) -> bool:
        v0 = Point(self.P3.X - self.P1.X, self.P3.Y - self.P1.Y, 0)
        v1 = Point(self.P2.X - self.P1.X, self.P2.Y - self.P1.Y, 0)
        v2 = Point(x - self.P1.X, y - self.P1.Y, 0)
        dot00 = v0.X ** 2 + v0.Y ** 2
        dot01 = v0.X * v1.X + v0.Y * v1.Y
        dot02 = v0.X * v2.X + v0.Y * v2.Y
        dot11 = v1.X ** 2 + v1.Y ** 2
        dot12 = v1.X * v2.X + v1.Y * v2.Y
        inv_denom = 1.0 / (dot00 * dot11 - dot01 ** 2)
        u = (dot11 * dot02 - dot01 * dot12) * inv_denom
        v = (dot00 * dot12 - dot01 * dot02) * inv_denom
        return (u > 0 or isclose(u, 0)) and (v > 0 or isclose(v, 0)) and (u + v < 1 or isclose(u + v, 1))

    def interpolate(self, x: float, y: float) -> float:
        ab = Point(self.P2.X - self.P1.X, self.P2.Y - self.P1.Y, self.P2.Z - self.P1.Z)
        ac = Point(self.P3.X - self.P1.X, self.P3.Y - self.P1.Y, self.P3.Z - self.P1.Z)
        X = ab.Y * ac.Z - ab.Z * ac.Y
        Y = ab.Z * ac.X - ab.X * ac.Z
        Z = ab.X * ac.Y - ab.Y * ac.X
        W = self.P1.X * X + self.P1.Y * Y + self.P1.Z * Z
        return (W - X * x - Y * y) / Z


class DelaunayMap:
    def __init__(self, *_points: Iterable[Tuple[float, float, float]]):
        points = sorted(list(Point(p[0], p[1], p[2]) for p in _points))
        if len(points) != len({(p.X, p.Y) for p in points}):
            raise Exception("(X,Y) coordinates of the points must be unique!")
        if len(points) < 3:
            raise Exception("There must be at least 3 points!")

        min_x, min_y = min(p.X for p in points), min(p.Y for p in points)
        max_x, max_y = max(p.X for p in points), max(p.Y for p in points)
        supertriangle = Triangle(Point(min_x - 1, min_y - 1, 0), Point(max_x * 3, min_y - 1, 0),
                                 Point(min_x - 1, max_y * 3, 0))
        triangles = [supertriangle]
        point_count = 0
        for p in points:
            point_count += 1
            container_triangles = [tr for tr in triangles if tr.contains_in_circumcircle(p)]
            edges = []  # type: List[Edge]
            for tr in container_triangles:
                triangles.remove(tr)
                edges += [Edge(tr.P1, tr.P2), Edge(tr.P2, tr.P3), Edge(tr.P3, tr.P1)]
            unique_edges = (edge for edge, count in Counter(edges).items() if count == 1)
            new_triangles = (Triangle(e.P1, e.P2, p) for e in unique_edges)
            triangles += new_triangles

            if len(triangles) != 2 * (point_count + 3) - 2 - 3:
                raise Exception(
                    "Triangulation invariant violated! Triangle count and point count doesn't fit together.")

        self._triangles = [tr for tr in triangles if not tr.shares_vertex(supertriangle)]

    def __getitem__(self, xy: Tuple[float, float]) -> float:
        x, y = xy
        for tr in self._triangles:  # type: Triangle
            if tr.is_inside(x, y):
                return tr.interpolate(x, y)


if __name__ == '__main__':
    with open("map.csv") as f:
        points = [tuple(float(x) for x in line.strip().split(",")) for line in f.readlines()]

    start = time.clock()
    d = DelaunayMap(*points)

    for x, y, z in points:
        assert z == d[x, y]
    print(time.clock() - start)
