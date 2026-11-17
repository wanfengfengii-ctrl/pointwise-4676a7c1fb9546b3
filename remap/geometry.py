"""Exact convex tetrahedron/closed-box intersection in global coordinates."""
from fractions import Fraction
from itertools import combinations
from .errors import ValidationError
from .models import Moments, ZERO, vector


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), Fraction(0))


def cross(a, b):
    return (a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0])


def determinant(a, b, c):
    return dot(a, cross(b, c))


def tetrahedron_moments(vertices):
    vertices = tuple(vector(p) for p in vertices)
    if len(vertices) != 4:
        raise ValidationError("expected four vertices")
    p, q, r, s = vertices
    volume = abs(determinant(sub(q, p), sub(r, p), sub(s, p))) / 6
    if not volume:
        raise ValidationError("degenerate source tetrahedron")
    return Moments(volume, tuple(volume * sum(p[i] for p in vertices) / 4
                                 for i in range(3)))


def _planes(vertices, lower, upper):
    planes = []
    for omitted in range(4):
        face = [vertices[i] for i in range(4) if i != omitted]
        p, q, r = face
        normal = cross(sub(q, p), sub(r, p))
        if dot(normal, sub(vertices[omitted], p)) > 0:
            normal = tuple(-x for x in normal)
        planes.append((normal, dot(normal, p)))
    for axis in range(3):
        n = tuple(Fraction(int(i == axis)) for i in range(3))
        planes.append((n, upper[axis]))
        planes.append((tuple(-x for x in n), -lower[axis]))
    return planes


def _meet(a, b, c):
    n, d = a
    m, e = b
    k, f = c
    det = determinant(n, m, k)
    if not det:
        return None
    u, v, w = cross(m, k), cross(k, n), cross(n, m)
    return tuple((d*u[i] + e*v[i] + f*w[i]) / det for i in range(3))


def _face_hull(points, normal):
    # Dropping a coordinate with nonzero normal is an injective plane projection.
    drop = next(i for i in range(3) if normal[i])
    axes = [i for i in range(3) if i != drop]
    points = sorted(points, key=lambda p: (p[axes[0]], p[axes[1]]))

    def turn(a, b, c):
        x, y = axes
        return (b[x]-a[x])*(c[y]-a[y]) - (b[y]-a[y])*(c[x]-a[x])

    def chain(seq):
        out = []
        for point in seq:
            while len(out) >= 2 and turn(out[-2], out[-1], point) <= 0:
                out.pop()
            out.append(point)
        return out

    return chain(points)[:-1] + chain(reversed(points))[:-1]


def intersection_moments(vertices, lower, upper):
    vertices = tuple(vector(p) for p in vertices)
    tetrahedron_moments(vertices)
    lower, upper = vector(lower), vector(upper)
    if any(a >= b for a, b in zip(lower, upper)):
        raise ValidationError("invalid target bounds")
    planes = _planes(vertices, lower, upper)
    points = set()
    for triple in combinations(planes, 3):
        p = _meet(*triple)
        if p is not None and all(dot(n, p) <= d for n, d in planes):
            points.add(p)
    points = sorted(points)
    if len(points) < 4:
        return ZERO
    origin = points[0]
    if not any(determinant(sub(a, origin), sub(b, origin), sub(c, origin))
               for a, b, c in combinations(points[1:], 3)):
        return ZERO
    center = tuple(sum(p[i] for p in points) / len(points) for i in range(3))
    volume = Fraction(0)
    first = [Fraction(0)] * 3
    seen = set()
    for normal, bound in planes:
        face = frozenset(p for p in points if dot(normal, p) == bound)
        if len(face) < 3 or face in seen:
            continue
        seen.add(face)
        hull = _face_hull(face, normal)
        if len(hull) < 3:
            continue
        a = hull[0]
        for i in range(1, len(hull)-1):
            b, c = hull[i:i+2]
            part = abs(determinant(sub(a, center), sub(b, center),
                                   sub(c, center))) / 6
            volume += part
            for axis in range(3):
                first[axis] += part * (center[axis]+a[axis]+b[axis]+c[axis]) / 4
    return Moments(volume, tuple(first))
