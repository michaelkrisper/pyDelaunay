from collections import Counter
from typing import List, Tuple

class Point:
    def __init__(self, x: float, y: float, z: float):
        self.Z = z
        self.X = x
        self.Y = y

    def __sub__(self, other: 'Point') -> 'Point':
        return Point(self.X - other.X, self.Y - other.Y, self.Z - other.Z)


class Edge:
    def __init__(self, p1: Point, p2: Point):
        self.P2 = p2
        self.P1 = p1


class Triangle:
    def __init__(self, p1: Point, p2: Point, p3: Point):
        self.P3 = p3
        self.P2 = p2
        self.P1 = p1

    def contains_in_circumcircle(self, p: Point) -> bool:
        p0 = self.P1 - p
        p1 = self.P2 - p
        p2 = self.P3 - p
        p0Square = p0.X ** 2 + p0.Y ** 2
        p1Square = p1.X ** 2 + p1.Y ** 2
        p2Square = p2.X ** 2 + p2.Y ** 2
        det01 = p0.X * p1.Y - p1.X * p0.Y
        det12 = p1.X * p2.Y - p2.X * p1.Y
        det20 = p2.X * p0.Y - p0.X * p2.Y
        return p0Square * det12 + p1Square * det20 + p2Square * det01 > 0

    def shares_vertex(self, tr: 'Triangle') -> bool:
        return any({self.P1, self.P2, self.P3} & {tr.P1, tr.P2, tr.P3})

    def is_inside(self, x: float, y: float) -> bool:
        if (self.P1.Y < y and self.P2.Y < y and self.P3.Y < y) or (
                            self.P1.X < x and self.P2.X < x and self.P3.X < x) or (
                            self.P1.X > x and self.P2.X > x and self.P3.X > x) or (
                            self.P1.Y > y and self.P2.Y > y and self.P3.Y > y):
            return False

        v0 = self.P3 - self.P1
        v1 = self.P2 - self.P1
        v2 = Point(x, y, 0) - self.P1

        dot00 = v0.X ** 2 + v0.Y ** 2
        dot01 = v0.X * v1.X + v0.Y * v1.Y
        dot02 = v0.X * v2.X + v0.Y * v2.Y
        dot11 = v1.X * v1.X + v1.Y * v1.Y
        dot12 = v1.X * v2.X + v1.Y * v2.Y

        inv_denom = dot00 * dot11 - dot01 ** 2
        u = (dot11 * dot02 - dot01 * dot12) / inv_denom
        v = (dot00 * dot12 - dot01 * dot02) / inv_denom

        return u >= 0 and v >= 0 and u + v <= 1


class Plane:
    def __init__(self, tr: Triangle):
        ab = tr.P2 - tr.P1
        ac = tr.P3 - tr.P1

        self.X = ab.Y * ac.Z - ab.Z * ac.Y
        self.Y = ab.Z * ac.X - ab.X * ac.Z
        self.Z = ab.X * ac.Y - ab.Y * ac.X
        self.W = tr.P1.X * self.X + tr.P1.Y * self.Y + tr.P1.Z * self.Z

    def __repr__(self) -> str:
        return "Plane({X},{Y},{Z},{W})".format_map(vars(self))


class DelaunayMap:
    def __init__(self, points: List[Point]):
        assert len(points) > 3

        min_x, min_y = min(p.X for p in points), min(p.Y for p in points)
        max_x, max_y = max(p.X for p in points), max(p.Y for p in points)

        supertriangle = Triangle(Point(min_x - 1, min_y - 1, 0), Point(max_x * 3, 0, 0), Point(min_x - 1, max_y * 3, 0))
        triangles = [supertriangle]

        for p in points:
            container_triangles = [tr for tr in triangles if tr.contains_in_circumcircle(p)]
            for tr in container_triangles:
                triangles.remove(tr)

            edges = [Edge(x, y)
                     for tr in container_triangles
                     for x, y in ((tr.P1, tr.P2), (tr.P2, tr.P3), (tr.P3, tr.P1))]
            unique_edges = [edge for edge, count in Counter(edges).items() if count == 1]
            triangles += [Triangle(e.P1, e.P2, p) for e in unique_edges]
        self.triangles = [tr for tr in triangles
                          if not tr.shares_vertex(supertriangle)]

    def __getitem__(self, xy: Tuple[float, float]) -> float:
        x, y = xy
        for tr in self.triangles:  # type: Triangle
            if tr.is_inside(x, y):
                pl = Plane(tr)
                return (pl.W - pl.X * x - pl.Y * y) / pl.Z

