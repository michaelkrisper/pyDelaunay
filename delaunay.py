from collections import Counter
from typing import List, Tuple, Dict
from math import isclose


class Point:
    def __init__(self, x: float, y: float, z: float):
        self.Z = z
        self.X = x
        self.Y = y


class Triangle:
    def __init__(self, p1: Point, p2: Point, p3: Point):
        self.P3 = p3
        self.P2 = p2
        self.P1 = p1

    def contains_in_circumcircle(self, p: Point):
        p0X = self.P1.X - p.X
        p0Y = self.P1.Y - p.Y
        p1X = self.P2.X - p.X
        p1Y = self.P2.Y - p.Y
        p2X = self.P3.X - p.X
        p2Y = self.P3.Y - p.Y
        p0Square = p0X * p0X + p0Y * p0Y
        p1Square = p1X * p1X + p1Y * p1Y
        p2Square = p2X * p2X + p2Y * p2Y
        det01 = p0X * p1Y - p1X * p0Y
        det12 = p1X * p2Y - p2X * p1Y
        det20 = p2X * p0Y - p0X * p2Y
        result = p0Square * det12 + p1Square * det20 + p2Square * det01
        return result > 0

    def shares_vertex(self, tr: 'Triangle'):
        return any({self.P1, self.P2, self.P3} & {tr.P1, tr.P2, tr.P3})

    def is_inside(self, x, y):
        if (self.P1.Y < y and self.P2.Y < y and self.P3.Y < y) or (
                            self.P1.X < x and self.P2.X < x and self.P3.X < x) or (
                            self.P1.X > x and self.P2.X > x and self.P3.X > x) or (
                            self.P1.Y > y and self.P2.Y > y and self.P3.Y > y):
            return False

        v0X = self.P3.X - self.P1.X
        v0Y = self.P3.Y - self.P1.Y
        v1X = self.P2.X - self.P1.X
        v1Y = self.P2.Y - self.P1.Y
        v2X = x - self.P1.X
        v2Y = y - self.P1.Y

        dot00 = v0X * v0X + v0Y * v0Y
        dot01 = v0X * v1X + v0Y * v1Y
        dot02 = v0X * v2X + v0Y * v2Y
        dot11 = v1X * v1X + v1Y * v1Y
        dot12 = v1X * v2X + v1Y * v2Y

        inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01)
        u = (dot11 * dot02 - dot01 * dot12) * inv_denom
        v = (dot00 * dot12 - dot01 * dot02) * inv_denom

        return u >= 0 and v >= 0 and u + v <= 1


class Edge:
    def __init__(self, p1: Point, p2: Point):
        self.P2 = p2
        self.P1 = p1


class Plane:
    def __init__(self, tr: Triangle):
        abX = tr.P2.X - tr.P1.X
        abY = tr.P2.Y - tr.P1.Y
        abZ = tr.P2.Z - tr.P1.Z

        acX = tr.P3.X - tr.P1.X
        acY = tr.P3.Y - tr.P1.Y
        acZ = tr.P3.Z - tr.P1.Z

        self.X = abY * acZ - abZ * acY
        self.Y = abZ * acX - abX * acZ
        self.Z = abX * acY - abY * acX
        self.W = tr.P1.X * self.X + tr.P1.Y * self.Y + tr.P1.Z * self.Z


class DelaunayMap:
    def __init__(self, points: List[Point]):
        assert len(points) > 3

        supertriangle = Triangle(Point(-1000, -1000, 0),
                                 Point(1000, 1000, 0),
                                 Point(1000, 0, 0))
        triangles = [supertriangle]

        for p in points:
            container_triangles = [tr for tr in triangles if tr.contains_in_circumcircle(p)]
            for tr in container_triangles:
                triangles.remove(tr)

            edges = [Edge(x, y)
                     for tr in container_triangles
                     for x, y in ((tr.P1, tr.P2), (tr.P2, tr.P3), (tr.P3, tr.P1))]
            c = Counter(edges)
            unique_edges = []
            for count, e in c.items():
                if count == 1:
                    unique_edges.append(e)

            new_triangles = [Triangle(e.P1, e.P2, p) for e in unique_edges]

            triangles += new_triangles
        self.triangles = [tr for tr in triangles
                          if not tr.shares_vertex(supertriangle)]

    def __getitem__(self, xy: Tuple[float, float]):
        x, y = xy
        for tr in self.triangles:  # type: Triangle
            if tr.is_inside(x, y):
                pl = Plane(tr)
                return x, y, (pl.W - pl.X * x - pl.Y * y) / pl.Z


if __name__ == '__main__':
    points = [Point(0, 0, 0),
              Point(1, 0, 1),
              Point(0, 1, 1),
              Point(1, 1, 1)]
    d = DelaunayMap(points)

    print(d[0.5, 0.5])
    assert d[0.5, 0.5] == 0.5

    print("finished")
