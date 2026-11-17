"""Exact planar-polygon tools used by the directed face evidence.

All coordinates and scalars are Fractions.  A polygon is a list of distinct
three-dimensional points in a single plane, in boundary order, without a
repeated first point.  No square roots are ever needed: the vector area of a
triangle, ``cross(b-a, c-a)/2``, already equals ``integral n dS`` because the
unit normal carries the same irrational factor as the area element.
"""
from fractions import Fraction

from .geometry import cross, dot, sub

ZERO3 = (Fraction(0), Fraction(0), Fraction(0))
MONOMIALS = ("1", "x", "y", "z", "xx", "xy", "xz", "yy", "yz", "zz")


def polygon_area_vector(poly):
    """Return integral of the outward unit normal over the polygon."""
    total = [Fraction(0), Fraction(0), Fraction(0)]
    if len(poly) < 3:
        return tuple(total)
    a = poly[0]
    for i in range(1, len(poly) - 1):
        c = cross(sub(poly[i], a), sub(poly[i + 1], a))
        for k in range(3):
            total[k] += c[k] / 2
    return tuple(total)


def _average_monomials(points):
    s = [sum((p[k] for p in points), Fraction(0)) for k in range(3)]
    return {
        "1": Fraction(1),
        "x": s[0] / 3, "y": s[1] / 3, "z": s[2] / 3,
        "xx": (sum((p[0] * p[0] for p in points), Fraction(0)) + s[0] * s[0]) / 12,
        "xy": (sum((p[0] * p[1] for p in points), Fraction(0)) + s[0] * s[1]) / 12,
        "xz": (sum((p[0] * p[2] for p in points), Fraction(0)) + s[0] * s[2]) / 12,
        "yy": (sum((p[1] * p[1] for p in points), Fraction(0)) + s[1] * s[1]) / 12,
        "yz": (sum((p[1] * p[2] for p in points), Fraction(0)) + s[1] * s[2]) / 12,
        "zz": (sum((p[2] * p[2] for p in points), Fraction(0)) + s[2] * s[2]) / 12,
    }


def polygon_moments(poly):
    """Exact ``integral m * n_i dS`` for the ten published monomials.

    A fan triangulation from the first vertex is used; the result is
    independent of the chosen fan for a planar polygon.
    """
    out = {m: [Fraction(0), Fraction(0), Fraction(0)] for m in MONOMIALS}
    if len(poly) < 3:
        return {m: tuple(v) for m, v in out.items()}
    a = poly[0]
    for i in range(1, len(poly) - 1):
        tri = (a, poly[i], poly[i + 1])
        area_vec = cross(sub(tri[1], tri[0]), sub(tri[2], tri[0]))
        if area_vec == ZERO3:
            continue
        avg = _average_monomials(tri)
        for m in MONOMIALS:
            weight = avg[m] / 2
            for k in range(3):
                out[m][k] += weight * area_vec[k]
    return {m: tuple(v) for m, v in out.items()}


def clip_polygon(poly, normal, bound):
    """Keep the part of a 3D planar polygon with ``dot(normal, p) <= bound``."""
    if not poly:
        return []
    out = []
    previous = poly[-1]
    previous_side = dot(normal, previous) - bound
    for point in poly:
        side = dot(normal, point) - bound
        if side <= 0:
            if previous_side > 0:
                out.append(_intersect(previous, point, previous_side, side))
            out.append(point)
        elif previous_side <= 0:
            out.append(_intersect(previous, point, previous_side, side))
        previous, previous_side = point, side
    return _dedupe(out)


def _intersect(p, q, sp, sq):
    # sp > 0 >= sq or sp <= 0 < sq; solve p + t (q-p) on the plane.
    t = sp / (sp - sq)
    return tuple(p[k] + t * (q[k] - p[k]) for k in range(3))


def _dedupe(points):
    out = []
    for p in points:
        if not out or p != out[-1]:
            out.append(p)
    while len(out) > 1 and out[0] == out[-1]:
        out.pop()
    return out


def clean_polygon(poly):
    """Remove consecutive duplicates and the repeated first point."""
    return _dedupe(list(poly))


# ---------------------------------------------------------------- 2D helpers

def _axes(axis):
    return tuple(i for i in range(3) if i != axis)


def _to2(poly, axis):
    a, b = _axes(axis)
    return [(p[a], p[b]) for p in poly]


def _unproject(points, axis, value):
    a, b = _axes(axis)
    out = []
    for u, v in points:
        p = [None, None, None]
        p[axis] = value
        p[a], p[b] = u, v
        out.append(tuple(p))
    return out


def _cross2(o, a, b):
    return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])


def _area2(poly):
    total = Fraction(0)
    for i in range(len(poly)):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % len(poly)]
        total += x1*y2 - x2*y1
    return total / 2


def _ccw(poly):
    if len(poly) >= 3 and _area2(poly) < 0:
        return list(reversed(poly))
    return list(poly)


def _clip2(subject, a, b):
    out = []
    if not subject:
        return out
    prev = subject[-1]
    prev_in = _cross2(a, b, prev) >= 0
    for cur in subject:
        cur_in = _cross2(a, b, cur) >= 0
        if cur_in:
            if not prev_in:
                out.append(_cross2d(prev, cur, a, b))
            out.append(cur)
        elif prev_in:
            out.append(_cross2d(prev, cur, a, b))
        prev, prev_in = cur, cur_in
    cleaned = []
    for p in out:
        if not cleaned or p != cleaned[-1]:
            cleaned.append(p)
    while len(cleaned) > 1 and cleaned[0] == cleaned[-1]:
        cleaned.pop()
    return cleaned


def _cross2d(p, q, a, b):
    num = _cross2(a, b, p)
    denominator = _cross2(a, b, p) - _cross2(a, b, q)
    t = num / denominator
    return (p[0] + t*(q[0]-p[0]), p[1] + t*(q[1]-p[1]))


def intersect_axis(poly_a, poly_b, axis, value):
    """Convex intersection of two coplanar polygons on ``p[axis] == value``."""
    if len(poly_a) < 3 or len(poly_b) < 3:
        return []
    subject, clip = _ccw(_to2(poly_a, axis)), _ccw(_to2(poly_b, axis))
    for i in range(len(clip)):
        subject = _clip2(subject, clip[i], clip[(i + 1) % len(clip)])
        if not subject:
            return []
    if _area2(subject) == 0:
        return []
    return _unproject(subject, axis, value)


def split_axis_polygon(poly, axis, value, line_a, line_b):
    """Split a coplanar polygon by the infinite line through two 3D points."""
    points = _to2(poly, axis)
    a2, b2 = (line_a[_axes(axis)[0]], line_a[_axes(axis)[1]]), (
        line_b[_axes(axis)[0]], line_b[_axes(axis)[1]])
    left = _unproject(_clip2(points, a2, b2), axis, value)
    right = _unproject(_clip2(points, b2, a2), axis, value)
    return clean_polygon(left), clean_polygon(right)


def decompose_outside(poly, axis, value, holes):
    """Cells of ``poly`` outside every hole polygon (holes have disjoint interiors).

    Holes are subtracted one at a time so the cell count stays linear in the
    number of holes instead of exponential in the number of hole edges.
    """
    cells = [poly]
    for hole in holes:
        next_cells = []
        for cell in cells:
            next_cells.extend(_subtract_convex(cell, hole, axis, value))
        cells = next_cells
    return cells


def _subtract_convex(cell, hole, axis, value):
    overlap = intersect_axis(cell, hole, axis, value)
    if not overlap or polygon_area_vector(overlap) == ZERO3:
        return [cell]
    pieces = [cell]
    for i in range(len(hole)):
        p, q = hole[i], hole[(i + 1) % len(hole)]
        if p == q:
            continue
        split = []
        for piece in pieces:
            left, right = split_axis_polygon(piece, axis, value, p, q)
            for part in (left, right):
                if len(part) >= 3 and polygon_area_vector(part) != ZERO3:
                    split.append(part)
        pieces = split
        if not pieces:
            return []
    result = []
    for piece in pieces:
        centroid = tuple(sum((p[k] for p in piece), Fraction(0)) / len(piece)
                         for k in range(3))
        if not _point_in_closed_convex(centroid, hole, axis, value):
            result.append(piece)
    return result


def _point_in_closed_convex(point, poly, axis, value):
    pts = _ccw(_to2(poly, axis))
    p2 = _to2([point], axis)[0]
    for i in range(len(pts)):
        a, b = pts[i], pts[(i + 1) % len(pts)]
        if _cross2(a, b, p2) < 0:
            return False
    return True
