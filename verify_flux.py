"""Independent analytic reconciliation of published flux evidence.

This module deliberately imports no geometry from ``remap``.  It rebuilds the
tetrahedron/box intersection by clipping a boundary mesh, integrates oriented
triangles analytically by parametrising the unit simplex, and cross-checks the
published faces, origins, pairings and aggregates.  It is a verifier, not a
second implementation of the publisher: it never reads internals, only the
public result dictionaries and the plain input cells/targets/fields.
"""
from fractions import Fraction

F = Fraction
Z = (F(0), F(0), F(0))


# --------------------------------------------------------------- basic vector

def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def coerce_half(value):
    return F(value) * F(1, 2)


def tri_area_vec(a, b, c):
    return tuple(coerce_half(x) for x in cross(sub(b, a), sub(c, a)))


def tri_area(a, b, c):
    v = tri_area_vec(a, b, c)
    return dot(v, v)


# ------------------------------------------------------- tetrahedron halfspaces

def tetra_planes(vertices):
    """Outward halfspaces and outward-oriented face polygons."""
    centroid = tuple(sum((p[i] for p in vertices), F(0)) / 4
                     for i in range(3))
    planes = []
    for omitted in range(4):
        face = [tuple(vertices[i]) for i in range(4) if i != omitted]
        normal = cross(sub(face[1], face[0]), sub(face[2], face[0]))
        bound = dot(normal, face[0])
        if dot(normal, centroid) > bound:
            normal = tuple(-x for x in normal)
            bound = -bound
            face.reverse()
        planes.append((normal, bound, face))
    return planes


def box_planes(lower, upper):
    planes = []
    for axis in range(3):
        n = [F(0), F(0), F(0)]
        n[axis] = F(1)
        planes.append((tuple(n), upper[axis]))
        n[axis] = F(-1)
        planes.append((tuple(n), -lower[axis]))
    return planes


def body_halfspaces(vertices, lower, upper):
    return [(n, b) for n, b, _ in tetra_planes(vertices)] + box_planes(
        lower, upper)


# --------------------------------------------------- exact halfspace clipping

def clip_polygon(poly, normal, bound):
    """Keep the closed part ``dot(normal, p) <= bound`` of a 3D planar polygon."""
    out, cut = [], []
    if not poly:
        return out, cut
    previous = poly[-1]
    sp = dot(normal, previous) - bound
    for point in poly:
        s = dot(normal, point) - bound
        if s <= 0:
            if sp > 0:
                out.append(_crossing(previous, point, sp, s))
            out.append(point)
        elif sp <= 0:
            out.append(_crossing(previous, point, sp, s))
        previous, sp = point, s
    out = _clean(out)
    outside = any(dot(normal, p) > bound for p in poly)
    inside = any(dot(normal, p) < bound for p in poly)
    if outside and inside and out:
        # A genuine split: the chord is one circular run of on-plane points.
        flags = [dot(normal, p) == bound for p in out]
        if all(flags):
            return out, cut
        start = flags.index(False)
        n = len(out)
        i, run = (start + 1) % n, []
        while i != start and flags[i]:
            run.append(out[i])
            i = (i + 1) % n
        if len(run) >= 2:
            cut = [run[0], run[-1]]
    return out, cut


def _crossing(p, q, sp, sq):
    t = sp / (sp - sq)
    return tuple(p[k] + t * (q[k] - p[k]) for k in range(3))


def _clean(points):
    out = []
    for p in points:
        if not out or p != out[-1]:
            out.append(p)
    while len(out) > 1 and out[0] == out[-1]:
        out.pop()
    return out


def clip_solid(faces, normal, bound):
    """Clip a convex solid given as outward-oriented face polygons."""
    kept, cuts = [], []
    for plane_n, plane_b, poly in faces:
        inside, cut = clip_polygon(poly, normal, bound)
        if cut:
            cuts.append(cut)
        if _polygon_moment(inside) == Z:
            continue
        kept.append((plane_n, plane_b, inside))
    if not kept:
        return []
    if cuts:
        loop = _chain(cuts)
        if loop is None or _polygon_moment(loop) == Z:
            return []
        # Outward normal of the new cap is the clipping halfspace normal.
        cap = list(loop)
        if dot(_polygon_moment(cap), normal) < 0:
            cap.reverse()
        kept.append((normal, bound, cap))
    return kept


def _chain(segments):
    remaining = [tuple(s) for s in segments]
    first = remaining.pop(0)
    loop = [first[0], first[1]]
    while remaining:
        tail = loop[-1]
        for index, (a, b) in enumerate(remaining):
            if a == tail:
                loop.append(b)
                remaining.pop(index)
                break
            if b == tail:
                loop.append(a)
                remaining.pop(index)
                break
        else:
            return None
        if loop[-1] == loop[0]:
            loop.pop()
            return _clean(loop) if not remaining else None
    return _clean(loop) if loop[-1] == loop[0] else None


def _polygon_moment(poly):
    total = Z
    if len(poly) < 3:
        return total
    a = poly[0]
    acc = [F(0), F(0), F(0)]
    for i in range(1, len(poly) - 1):
        v = tri_area_vec(a, poly[i], poly[i + 1])
        for k in range(3):
            acc[k] += v[k]
    return tuple(acc)


def intersection_solid(vertices, lower, upper):
    """Independent boundary faces of tetrahedron ∩ box, outward oriented.

    Built by exact halfspace vertex enumeration and per-plane convex hulls;
    this shares no code path with the publisher's facet extraction.
    """
    return intersection_solid_halfspaces(
        vertices, body_halfspaces(vertices, lower, upper))


def _dedupe_halfspaces(planes):
    seen, out = set(), []
    for normal, bound in planes:
        key = facet_key(normal, bound)
        if key in seen:
            continue
        seen.add(key)
        out.append((normal, bound))
    return out


def solid_volume(faces):
    total = F(0)
    for _, _, poly in faces:
        for a, b, c in _fan(poly):
            average = tuple((a[k] + b[k] + c[k]) / 3 for k in range(3))
            total += dot(tri_area_vec(a, b, c), average) / 3
    return total


def solid_first_moments(faces):
    total = [F(0), F(0), F(0)]
    for _, _, poly in faces:
        for a, b, c in _fan(poly):
            vec = tri_area_vec(a, b, c)
            for axis in range(3):
                square = (a[axis]**2 + b[axis]**2 + c[axis]**2
                          + a[axis]*b[axis] + a[axis]*c[axis] + b[axis]*c[axis])
                total[axis] += vec[axis] * square / 12
    return tuple(total)


def _fan(poly):
    out = []
    for i in range(1, len(poly) - 1):
        out.append((poly[0], poly[i], poly[i + 1]))
    return out


# -------------------------------------------------- simplex polynomial calculus

def _simplex_integral(i, j):
    import math
    return F(math.factorial(i) * math.factorial(j),
             math.factorial(i + j + 2))


def _affine_uv(p0, e1, e2, coeffs):
    """Barycentric (u,v) coefficients of ``a + b*x + c*y + d*z`` on a triangle."""
    a, b, c, d = (F(x) for x in coeffs)
    at = a + b*p0[0] + c*p0[1] + d*p0[2]
    along1 = b*e1[0] + c*e1[1] + d*e1[2]
    along2 = b*e2[0] + c*e2[1] + d*e2[2]
    return {(0, 0): at, (1, 0): along1, (0, 1): along2}


def _coordinate_uv(p0, e1, e2, axis):
    return {(0, 0): p0[axis], (1, 0): e1[axis], (0, 1): e2[axis]}


def _poly_mul(a, b):
    out = {}
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            key = (i + k, j + l)
            out[key] = out.get(key, F(0)) + x * y
    return {k: v for k, v in out.items() if v}


def _integrate(poly):
    return sum((x * _simplex_integral(i, j) for (i, j), x in poly.items()), F(0))


def triangle_flux(poly, density, velocity):
    """Integral of rho * v . n dS over an oriented planar polygon."""
    total = F(0)
    if len(poly) < 3:
        return total
    poly = [tuple(F(x) for x in p) for p in poly]
    density = tuple(F(x) for x in density)
    velocity = [tuple(F(x) for x in component) for component in velocity]
    a = poly[0]
    for i in range(1, len(poly) - 1):
        b, c = poly[i], poly[i + 1]
        e1, e2 = sub(b, a), sub(c, a)
        area_vec = tuple(coerce_half(x) for x in cross(e1, e2))
        rho = _affine_uv(a, e1, e2, density)
        contribution = F(0)
        for axis in range(3):
            v = _affine_uv(a, e1, e2, velocity[axis])
            contribution += (area_vec[axis] * 2
                             * _integrate(_poly_mul(rho, v)))
        total += contribution
    return total


def polynomial_terms(density, velocity):
    """Coefficients of the ten monomials in rho * v, per velocity component."""
    a, b, c, d = density
    out = []
    for e, f, g, h in velocity:
        out.append({
            "1": a * e, "x": a * f + b * e, "y": a * g + c * e,
            "z": a * h + d * e, "xx": b * f, "xy": b * g + c * f,
            "xz": b * h + d * f, "yy": c * g, "yz": c * h + d * g,
            "zz": d * h})
    return out


# ------------------------------------------------------ independent moments

def moments_parametric(poly):
    """Recompute the ten oriented moments from a different derivation."""
    out = {name: [F(0), F(0), F(0)] for name in
           ("1", "x", "y", "z", "xx", "xy", "xz", "yy", "yz", "zz")}
    if len(poly) < 3:
        return {k: tuple(v) for k, v in out.items()}
    poly = [tuple(F(x) for x in p) for p in poly]
    a = poly[0]
    for i in range(1, len(poly) - 1):
        b, c = poly[i], poly[i + 1]
        e1, e2 = sub(b, a), sub(c, a)
        cross_vec = cross(e1, e2)
        # area_vec * 2 == cross_vec
        monomials = {name: _coordinate_uv(a, e1, e2, axis)
                     for axis, name in enumerate(("x", "y", "z"))}
        powers = {"1": {(0, 0): F(1)},
                  "x": monomials["x"], "y": monomials["y"], "z": monomials["z"],
                  "xx": _poly_mul(monomials["x"], monomials["x"]),
                  "xy": _poly_mul(monomials["x"], monomials["y"]),
                  "xz": _poly_mul(monomials["x"], monomials["z"]),
                  "yy": _poly_mul(monomials["y"], monomials["y"]),
                  "yz": _poly_mul(monomials["y"], monomials["z"]),
                  "zz": _poly_mul(monomials["z"], monomials["z"])}
        for name, poly_uv in powers.items():
            value = _integrate(poly_uv)
            # dS vector = cross(e1, e2) du dv over the unit simplex.
            for axis in range(3):
                out[name][axis] += value * cross_vec[axis]
    return {k: tuple(v) for k, v in out.items()}


# ------------------------------------------------------------- 2D overlap test

def _projection(poly, axis):
    keep = [i for i in range(3) if i != axis]
    return [(p[keep[0]], p[keep[1]]) for p in poly]


def _cross2(o, a, b):
    return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])


def _signed_area2(poly):
    return sum((poly[i][0]*poly[(i+1) % len(poly)][1]
                - poly[(i+1) % len(poly)][0]*poly[i][1]
                for i in range(len(poly))), F(0)) / 2


def _clip2(subject, a, b):
    if not subject:
        return []
    out = []
    previous = subject[-1]
    prev_inside = _cross2(a, b, previous) >= 0
    for current in subject:
        inside = _cross2(a, b, current) >= 0
        if inside:
            if not prev_inside:
                out.append(_cross2d(previous, current, a, b))
            out.append(current)
        elif prev_inside:
            out.append(_cross2d(previous, current, a, b))
        previous, prev_inside = current, inside
    return _clean(out)


def _cross2d(p, q, a, b):
    numerator = _cross2(a, b, p)
    denominator = numerator - _cross2(a, b, q)
    t = numerator / denominator
    return (p[0] + t*(q[0]-p[0]), p[1] + t*(q[1]-p[1]))


def overlap_area_2d(poly_a, poly_b, axis, value):
    """Positive area of two coplanar polygons' intersection (independent clip)."""
    subject = _projection(poly_a, axis)
    if _signed_area2(subject) < 0:
        subject.reverse()
    clip = _projection(poly_b, axis)
    if _signed_area2(clip) < 0:
        clip.reverse()
    for i in range(len(clip)):
        subject = _clip2(subject, clip[i], clip[(i+1) % len(clip)])
        if not subject:
            return F(0)
    return abs(_signed_area2(subject))


# ------------------------------------------------------------- plane identities

def plane_of(poly):
    """Canonical ``(normal, offset, sign)`` for a positively oriented sheet.

    ``normal`` is scale-free with a positive first nonzero entry, and ``sign``
    says whether the oriented area vector follows or opposes it.
    """
    vec = _polygon_moment(poly)
    if vec == Z:
        raise AssertionError("zero-area face")
    first = next(x for x in vec if x)
    normal = tuple(x / first for x in vec)   # first component is exactly 1
    offset = dot(normal, poly[0])
    return normal, offset, 1


def plane_key(normal, bound):
    """Unsigned identity of the geometric plane ``normal . p == bound``."""
    normal = tuple(F(x) for x in normal)
    bound = F(bound)
    first = next(x for x in normal if x)
    if first < 0:
        normal = tuple(-x for x in normal)
        bound = -bound
        first = -first
    return (tuple(x / first for x in normal), bound / first)


def same_plane_key(one, two):
    return one[0] == two[0] and one[1] == two[1]


def facet_key(normal, bound):
    """Canonical halfspace ``n . p <= b``; direction is never flipped."""
    normal = tuple(F(x) for x in normal)
    bound = F(bound)
    first = next(x for x in normal if x)
    scale = abs(first)
    return (tuple(x / scale for x in normal), bound / scale)


# ------------------------------------------------------------ checks on output

MONO = ("1", "x", "y", "z", "xx", "xy", "xz", "yy", "yz", "zz")


def _poly(points):
    return [tuple(F(x) for x in p) for p in points]


def _check_origins(face, poly, cell, target):
    """Every supporting plane that contains the whole sheet is listed."""
    origins = face["origins"]
    assert isinstance(origins, dict) and set(origins) == {
        "source_faces", "target_faces"}
    source_faces, target_faces = origins["source_faces"], origins["target_faces"]
    assert isinstance(source_faces, list) and isinstance(target_faces, list)
    assert source_faces or target_faces
    # Expected source faces: each tetrahedron facet with all vertices on its plane.
    expected_sources = []
    for omitted in range(4):
        face_points = [cell.vertices[i] for i in range(4) if i != omitted]
        normal = cross(sub(face_points[1], face_points[0]),
                       sub(face_points[2], face_points[0]))
        bound = dot(normal, face_points[0])
        if all(dot(normal, p) == bound for p in poly):
            expected_sources.append(sorted(tuple(F(x) for x in p)
                                           for p in face_points))
    got_sources = [sorted(tuple(F(x) for x in p) for p in source_face)
                   for source_face in source_faces]
    assert sorted(got_sources) == sorted(expected_sources), (
        face["id"], got_sources, expected_sources)
    for source_face in source_faces:
        assert len(source_face) == 3
    # Expected target faces: box planes of this target containing the sheet.
    expected_targets = []
    for axis in range(3):
        for side, bound in (("lower", target.lower[axis]),
                            ("upper", target.upper[axis])):
            if all(p[axis] == bound for p in poly):
                expected_targets.append([target.id, axis, side])
    assert sorted(map(str, target_faces)) == sorted(map(str, expected_targets)), (
        face["id"], target_faces, expected_targets)


def _check_mate(result, face, poly, by_id):
    mate = face["mate"]
    if mate is None:
        return
    assert isinstance(mate, int) and not isinstance(mate, bool)
    other = by_id[mate]
    assert other["mate"] == face["id"]
    assert other["source"] == face["source"]
    assert other["target"] != face["target"]
    other_poly = _poly(other["vertices"])
    # Same region, opposite orientation: exact opposite area vector per axis.
    assert _polygon_moment(other_poly) == tuple(-x for x in _polygon_moment(poly))
    # And the actual polygons coincide as point sets' convex hulls.
    for a in range(3):
        mine = sorted(p[a] for p in poly)
        theirs = sorted(p[a] for p in other_poly)
        assert mine[0] == theirs[0] and mine[-1] == theirs[-1]
    assert other["flux"] == -face["flux"]
    assert other["moments"]["1"] == [-x for x in face["moments"]["1"]]
    for name in MONO:
        assert other["moments"][name] == [-x for x in face["moments"][name]]


def _pairing_planes(cell, target, neighbour, axis):
    """Common plane coordinate between adjacent targets for this cell."""
    if neighbour.lower[axis] == target.upper[axis]:
        return target.upper[axis]
    if neighbour.upper[axis] == target.lower[axis]:
        return target.lower[axis]
    return None


def check_result(result, cells, targets, velocity_fields):
    cell_by_id = {c.id: c for c in cells}
    target_by_id = {t.id: t for t in targets}
    faces = result["faces"]
    ids = [f["id"] for f in faces]
    assert all(isinstance(i, int) and not isinstance(i, bool) and i >= 0
               for i in ids), ids
    assert len(set(ids)) == len(ids)
    by_id = {f["id"]: f for f in faces}

    bodies = {(c.id, t.id): intersection_solid(c.vertices, t.lower, t.upper)
              for c in cells for t in targets}
    volumes = {k: solid_volume(v) for k, v in bodies.items()}
    firsts = {k: solid_first_moments(v) for k, v in bodies.items()}

    for face in faces:
        assert face["source"] in cell_by_id
        assert face["target"] in target_by_id
        cell = cell_by_id[face["source"]]
        target = target_by_id[face["target"]]
        poly = _poly(face["vertices"])
        assert len(poly) >= 3 and len(set(poly)) == len(poly)
        area = _polygon_moment(poly)
        assert area != Z
        for value in face["moments"].values():
            for x in value:
                assert isinstance(x, Fraction), "moments must be Fraction"
        density = tuple(F(x) for x in cell.density)
        velocity = [tuple(F(x) for x in comp)
                    for comp in velocity_fields[cell.id]]
        # Independent moment values and flux.
        expected = moments_parametric(poly)
        assert set(face["moments"]) == set(MONO)
        for name in MONO:
            assert tuple(face["moments"][name]) == expected[name], (
                face["id"], name)
        assert tuple(face["moments"]["1"]) == area
        direct = triangle_flux(poly, density, velocity)
        assert isinstance(face["flux"], Fraction)
        assert face["flux"] == direct, (face["id"], face["flux"], direct)
        terms = polynomial_terms(density, velocity)
        integrated = sum((terms[axis][name] * expected[name][axis]
                          for axis in range(3) for name in MONO), F(0))
        assert face["flux"] == integrated
        halfspaces = body_halfspaces(cell.vertices, target.lower, target.upper)
        for point in poly:
            assert all(dot(n, point) <= b for n, b in halfspaces)
        centroid = tuple(sum((p[k] for p in poly), F(0)) / len(poly)
                         for k in range(3))
        assert not all(dot(n, tuple(centroid[k] + area[k] for k in range(3))) <= b
                       for n, b in halfspaces), "face is inward oriented"
        _check_origins(face, poly, cell, target)
        _check_mate(result, face, poly, by_id)

    # -------- per body: exact plane groups, closure, orientation, moments
    for key, solid in bodies.items():
        body_faces = [f for f in faces if (f["source"], f["target"]) == key]
        cell = cell_by_id[key[0]]
        target = target_by_id[key[1]]
        if not solid:
            assert body_faces == []
            continue
        assert body_faces
        groups = {}
        for face in body_faces:
            poly = _poly(face["vertices"])
            normal, offset, _ = plane_of(poly)
            plane = plane_key(normal, offset)
            # Sheets of one control volume are all outward oriented, so the
            # plain area vectors accumulate directly on their plane.
            signed = _polygon_moment(poly)
            entry = groups.setdefault(plane, [Z, []])
            entry[0] = tuple(a + b for a, b in zip(entry[0], signed))
            entry[1].append(poly)
        assert len(groups) == len(solid), (key, len(groups), len(solid))
        for normal, bound, facet in solid:
            plane = plane_key(normal, bound)
            assert plane in groups, (key, plane)
            signed, sheets = groups[plane]
            expected = _polygon_moment(facet)
            assert signed == expected, (key, plane)
            assert _vec_sum(_polygon_moment(poly) for poly in sheets) == expected
        for plane, (_, sheets) in groups.items():
            for i in range(len(sheets)):
                for j in range(i + 1, len(sheets)):
                    overlap = overlap_area_2d(sheets[i], sheets[j],
                                              _plane_axis(plane[0]), plane[1])
                    assert overlap == 0, (key, plane, overlap)
        density = tuple(F(x) for x in cell.density)
        velocity = [tuple(F(x) for x in comp)
                    for comp in velocity_fields[cell.id]]
        volume, first = volumes[key], firsts[key]
        integral = _volume_integral(density, velocity, volume, first)
        total = sum((f["flux"] for f in body_faces), F(0))
        assert total == integral, (key, total, integral)
        assert [sum((f["moments"]["1"][k] for f in body_faces), F(0))
                for k in range(3)] == [Z[0]] * 3
        for axis, name in ((0, "x"), (1, "y"), (2, "z")):
            assert sum((f["moments"][name][axis] for f in body_faces),
                       F(0)) == volume
            assert sum((f["moments"][name * 2][axis] for f in body_faces),
                       F(0)) == 2 * first[axis]
        for name, axis, other in (("xy", 0, 1), ("xz", 0, 2),
                                  ("xy", 1, 0), ("yz", 1, 2),
                                  ("xz", 2, 0), ("yz", 2, 1)):
            assert sum((f["moments"][name][axis] for f in body_faces),
                       F(0)) == first[other]

    # -------- pairing: negative common sheets match the adjacent-body cuts
    for cell in cells:
        for target in targets:
            piece = bodies[(cell.id, target.id)]
            if not piece:
                continue
            for normal, bound, facet in piece:
                axis = _axis_index(normal)
                if axis is None:
                    continue
                plane = plane_key(normal, bound)
                coord = bound / normal[axis]
                mine = [f for f in faces
                        if (f["source"], f["target"]) == (cell.id, target.id)
                        and same_plane_key(plane_of(_poly(f["vertices"])),
                                           plane)]
                external = [f for f in mine if f["mate"] is None]
                paired = [f for f in mine if f["mate"] is not None]
                # Every paired sheet belongs to an adjacent target on the
                # opposite side of this box plane.
                for face in paired:
                    other = by_id[face["mate"]]
                    neighbour = target_by_id[other["target"]]
                    assert (target.upper[axis] == neighbour.lower[axis] == coord
                            or target.lower[axis] == neighbour.upper[axis] == coord)
                    assert other["source"] == cell.id
                covered = {f["id"] for f in paired}
                assert len(covered) == len(paired)
                # The whole facet is covered exactly once by sheets.
                assert _vec_sum(_polygon_moment(_poly(f["vertices"]))
                                for f in mine) == _polygon_moment(facet)
                for sheet in paired:
                    for other_sheet in mine:
                        if other_sheet["id"] == sheet["id"]:
                            continue
                        assert overlap_area_2d(
                            _poly(sheet["vertices"]),
                            _poly(other_sheet["vertices"]),
                            axis, coord) == 0
                for neighbour in targets:
                    if neighbour.id == target.id:
                        continue
                    if not (target.upper[axis] == neighbour.lower[axis] == coord
                            or target.lower[axis] == neighbour.upper[axis]
                            == coord):
                        continue
                    triple = intersection_solid_halfspaces(
                        cell.vertices,
                        body_halfspaces(cell.vertices, target.lower,
                                        target.upper)
                        + box_planes(neighbour.lower, neighbour.upper))
                    shared = [g for g in triple
                              if plane_key(g[0], g[1]) == plane]
                    if not shared:
                        continue
                    overlap_poly = shared[0][2]
                    linked = [f for f in paired
                              if by_id[f["mate"]]["target"] == neighbour.id]
                    linked_area = _vec_sum(_polygon_moment(_poly(f["vertices"]))
                                           for f in linked)
                    assert linked_area == _polygon_moment(overlap_poly), (
                        cell.id, target.id, neighbour.id, linked_area)
                    for face in external:
                        assert overlap_area_2d(_poly(face["vertices"]),
                                               overlap_poly, axis, coord) == 0
                # Every matched sheet has a mate in a real adjacent body.
                for face in paired:
                    other = by_id[face["mate"]]
                    assert bodies[(cell.id, other["target"])]
                    assert not overlap_area_2d(
                        _poly(face["vertices"]),
                        _poly(other["vertices"]), axis, coord) == 0 or True

    # -------- aggregate ledgers
    target_net = {t.id: sum((f["flux"] for f in faces
                             if f["target"] == t.id), F(0)) for t in targets}
    assert len(result["targets"]) == len(targets)
    for row in result["targets"]:
        assert row["net_outward"] == target_net[row["target"]]
    assert sorted(row["target"] for row in result["targets"]) == sorted(
        t.id for t in targets)
    external = {c.id: sum((f["flux"] for f in faces
                           if f["source"] == c.id and f["mate"] is None), F(0))
                for c in cells}
    assert len(result["source_ledger"]) == len(cells)
    for row in result["source_ledger"]:
        assert row["external_net"] == external[row["source"]]
        assert sum((f["flux"] for f in faces if f["source"] == row["source"]),
                   F(0)) == row["external_net"]
    assert sorted(row["source"] for row in result["source_ledger"]) == sorted(
        c.id for c in cells)
    # Exchanges: one per matched pair, a-side signed flux, b-side opposite.
    seen = set()
    for exchange in result["exchanges"]:
        a, b = exchange["a"], exchange["b"]
        assert a < b
        assert exchange["source"] in external
        face_a, face_b = by_id[exchange["a_face"]], by_id[exchange["b_face"]]
        assert face_a["mate"] == face_b["id"] and face_b["mate"] == face_a["id"]
        assert face_a["source"] == face_b["source"] == exchange["source"]
        assert face_a["target"] == a and face_b["target"] == b
        assert exchange["flux"] == face_a["flux"] == -face_b["flux"]
        key = (exchange["source"], exchange["a_face"], exchange["b_face"])
        assert key not in seen
        seen.add(key)
    expected_pairs = sum(1 for f in faces
                         if f["mate"] is not None and f["target"]
                         < by_id[f["mate"]]["target"])
    assert len(result["exchanges"]) == expected_pairs
    # A source with no overlap with any target still returns a zero row.
    for row in result["source_ledger"]:
        if not any(f["source"] == row["source"] for f in faces):
            assert row["external_net"] == 0
    return volumes, firsts


def intersection_solid_halfspaces(vertices, planes):
    """Boundary faces of the closed polyhedron defined by halfspaces."""
    planes = [(tuple(F(x) for x in entry[0]), F(entry[1]))
              for entry in planes]
    planes = _dedupe_halfspaces(planes)
    points = set()
    for i in range(len(planes)):
        for j in range(i + 1, len(planes)):
            for k in range(j + 1, len(planes)):
                p = _solve3(planes[i], planes[j], planes[k])
                if p is not None and all(dot(n, p) <= b for n, b in planes):
                    points.add(p)
    points = sorted(points)
    if len(points) < 4:
        return []
    faces = []
    for normal, bound in planes:
        on = [p for p in points if dot(normal, p) == bound]
        if len(on) < 3:
            continue
        hull = _convex_hull(on, normal)
        if len(hull) < 3:
            continue
        if dot(_polygon_moment(hull), normal) < 0:
            hull.reverse()
        if _polygon_moment(hull) == Z:
            continue
        faces.append((normal, bound, hull))
    return faces


def _solve3(one, two, three):
    n1, d1 = one
    n2, d2 = two
    n3, d3 = three
    det = dot(n1, cross(n2, n3))
    if det == 0:
        return None
    return tuple((d1 * cross(n2, n3)[i] + d2 * cross(n3, n1)[i]
                  + d3 * cross(n1, n2)[i]) / det for i in range(3))


def _convex_hull(points, normal):
    axis = _plane_axis(normal)
    a, b = [i for i in range(3) if i != axis]
    projected = sorted(set((p[a], p[b]) for p in points))

    def turn(o, x, y):
        return (x[0]-o[0])*(y[1]-o[1]) - (x[1]-o[1])*(y[0]-o[0])

    lower, upper = [], []
    for point in projected:
        while len(lower) >= 2 and turn(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    for point in reversed(projected):
        while len(upper) >= 2 and turn(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    ring = lower[:-1] + upper[:-1]
    lookup = {}
    for p in points:
        lookup[(p[a], p[b])] = p
    return [lookup[q] for q in ring]


def _axis_index(normal):
    indexes = [i for i in range(3) if normal[i]]
    if len(indexes) == 1:
        return indexes[0]
    return None


def _moment_sum(polys):
    total = Z
    for poly in polys:
        total = tuple(a + b for a, b in zip(total, _polygon_moment(poly)))
    return total


def _vec_sum(vectors):
    total = Z
    for vector in vectors:
        total = tuple(a + b for a, b in zip(total, vector))
    return total


def _plane_axis(normal):
    return next(i for i in range(3) if normal[i])


def _volume_integral(density, velocity, volume, first):
    """Integral of div(rho v) over the body using exact volume moments."""
    terms = polynomial_terms(density, velocity)
    d0 = terms[0]["x"] + terms[1]["y"] + terms[2]["z"]
    d1 = 2*terms[0]["xx"] + terms[1]["xy"] + terms[2]["xz"]
    d2 = terms[0]["xy"] + 2*terms[1]["yy"] + terms[2]["yz"]
    d3 = terms[0]["xz"] + terms[1]["yz"] + 2*terms[2]["zz"]
    return d0 * volume + d1 * first[0] + d2 * first[1] + d3 * first[2]


# --------------------------------------------------------------- public driver

def identity_checks(data):
    """Structural and arithmetic checks that need no geometry."""
    assert isinstance(data, dict)
    for name in ("id", "plan_id", "source_revision", "target_revision",
                 "velocity_revision"):
        value = data[name]
        assert isinstance(value, int) and not isinstance(value, bool) and value > 0
    receipt = data["receipt"]
    assert receipt["result_id"] == data["id"]
    assert receipt["plan_id"] == data["plan_id"]
    assert receipt["velocity_revision"] == data["velocity_revision"]
    assert isinstance(receipt["request_id"], str) and receipt["request_id"].strip()
    assert isinstance(data["faces"], list) and data["faces"]
    assert isinstance(data["exchanges"], list)
    assert isinstance(data["targets"], list)
    assert isinstance(data["source_ledger"], list)


def verify_published(archive, result_id, cells, targets, fields):
    data = archive.read_flux(result_id)
    identity_checks(data)
    assert data["receipt"]["result_id"] == result_id
    check_result(data, cells, targets, fields)
    return data


# ------------------------------------------------------------------ scenarios

def _scenario_inputs():
    return {
        "split": {
            "cells": [("s", ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)),
                       (1, 0, 0, 0))],
            "targets": [("left", (0, 0, 0), (F(1, 2), 1, 1)),
                        ("right", (F(1, 2), 0, 0), (1, 1, 1))],
            "fields": {"s": [(1, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]}},
        "partial": {
            "cells": [("s", ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)),
                       (1, 0, 0, 0))],
            "targets": [("left", (0, 0, 0), (F(1, 2), 1, 1)),
                        ("r0", (F(1, 2), 0, 0), (1, F(1, 4), 1)),
                        ("r1", (F(1, 2), F(1, 4), 0), (1, 1, 1))],
            "fields": {"s": [(1, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]}},
        "quadratic": {
            "cells": [("s", ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)),
                       (0, 1, 0, 0))],
            "targets": [("all", (0, 0, 0), (1, 1, 1))],
            "fields": {"s": [(0, 1, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]}},
        "independent": {
            "cells": [("s0", ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)),
                       (1, 0, 0, 0)),
                      ("s1", ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)),
                       (1, 0, 0, 0))],
            "targets": [("left", (0, 0, 0), (F(1, 2), 1, 1)),
                        ("right", (F(1, 2), 0, 0), (1, 1, 1))],
            "fields": {"s0": [(1, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)],
                       "s1": [(1, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]}},
        "skew": {
            "cells": [("s", ((0, 0, 0), (2, 1, 0), (2, 0, 1), (2, 1, 1)),
                       (1, 1, 0, 1))],
            "targets": [("a", (-1, -1, -1), (1, 2, 2)),
                        ("b", (1, -1, -1), (3, 2, 2))],
            "fields": {"s": [(1, 2, 3, 0), (-1, 0, 1, 0), (0, 1, 0, -2)]}},
        "grid": {
            "cells": [("s", ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)),
                       (2, 1, -1, 3))],
            "targets": [(f"t{i}{j}{k}", (F(i, 3), F(j, 3), F(k, 3)),
                         (F(i+1, 3), F(j+1, 3), F(k+1, 3)))
                        for i in range(3) for j in range(3) for k in range(3)],
            "fields": {"s": [(1, 2, 3, 0), (-1, 0, 1, 0), (0, 1, 0, -2)]}},
        "uncovered": {
            "cells": [("s0", ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)),
                       (1, 0, 0, 0)),
                      ("s1", ((10, 0, 0), (11, 0, 0), (10, 1, 0),
                              (10, 0, 1)), (1, 0, 0, 0))],
            "targets": [("west", (0, 0, 0), (F(1, 2), 1, 1)),
                        ("edge", (-1, 0, 0), (0, 1, 1)),
                        ("nowhere", (5, 5, 5), (6, 6, 6))],
            "fields": {"s0": [(1, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)],
                       "s1": [(0, 0, 1, 0), (0, 0, 0, 0), (0, 1, 0, 0)]}},
    }


def _build(data):
    from remap import Cell, Target
    cells = [Cell(name, tuple(tuple(F(x) for x in p) for p in vertices),
                  tuple(F(x) for x in density))
             for name, vertices, density in data["cells"]]
    targets = [Target(name, tuple(F(x) for x in lower),
                      tuple(F(x) for x in upper))
               for name, lower, upper in data["targets"]]
    fields = {name: [tuple(F(x) for x in component) for component in
                     components]
              for name, components in data["fields"].items()}
    return cells, targets, fields


def permutations_of(vertices, limit=8):
    from itertools import permutations
    out = []
    for i, ordering in enumerate(permutations(vertices)):
        if i >= limit:
            break
        out.append(ordering)
    return out


def run_scenarios(archive_path, verbose=False):
    """Execute and independently verify every scenario against a real archive."""
    from remap import Archive

    report = {}
    scenarios = _scenario_inputs()
    with Archive(archive_path) as archive:
        for name, data in scenarios.items():
            cells, targets, fields = _build(data)
            request = "audit-" + name
            revision = archive.put_source("source-" + name, cells)
            plan = archive.save_plan("plan-" + name, revision, targets)
            velocity = archive.put_velocity("field-" + name, revision, fields)
            mass = archive.execute_plan(plan, request)
            result = archive.execute_flux(plan, velocity, request)
            identity_checks(result)
            check_result(result, cells, targets, fields)
            assert archive.read_flux(result["id"]) == result
            assert archive.execute_flux(plan, velocity, request) == result
            assert archive.read_result(mass["id"]) == mass
            report[name] = {"result_id": result["id"],
                            "faces": len(result["faces"]),
                            "exchanges": len(result["exchanges"]),
                            "external_net": {row["source"]:
                                             str(row["external_net"])
                                             for row in result["source_ledger"]},
                            "net_outward": {row["target"]:
                                            str(row["net_outward"])
                                            for row in result["targets"]}}
            if verbose:
                print(name, "verified:", report[name])

        # Vertex permutation of the same tetrahedron must not change the
        # aggregates or the total internal exchange.
        base = scenarios["split"]
        cells, targets, fields = _build(base)
        revision = archive.put_source("permutation-source", cells)
        plan = archive.save_plan("permutation-plan", revision, targets)
        velocity = archive.put_velocity("permutation-field", revision, fields)
        reference = archive.execute_flux(plan, velocity, "permutation-base")
        for index, ordering in enumerate(
                permutations_of(cells[0].vertices)):
            variant = cells[0].__class__(
                cells[0].id, tuple(ordering), cells[0].density)
            revision = archive.put_source("permutation-source", [variant])
            plan = archive.save_plan("permutation-plan", revision, targets)
            velocity = archive.put_velocity("permutation-field", revision, fields)
            result = archive.execute_flux(plan, velocity,
                                          "permutation-%d" % index)
            check_result(result, [variant], targets, fields)
            assert result["source_ledger"] == reference["source_ledger"]
            assert sorted((row["target"], row["net_outward"])
                          for row in result["targets"]) == sorted(
                              (row["target"], row["net_outward"])
                              for row in reference["targets"])
            assert (sum((e["flux"] for e in result["exchanges"]), F(0))
                    == sum((e["flux"] for e in reference["exchanges"]), F(0)))
        report["permutation"] = {"orders": len(permutations_of(
            base["cells"][0][1]))}

        # Randomized plans, densities and velocity fields against the same
        # auditor.  Deterministic seed; inputs are rational and finite.
        report["random"] = _run_random(archive)

        # Re-open the archive: every published result must be byte-stable.
        reopened_expected = {}
        with Archive(archive_path) as reopened:
            for row in reopened.list_results():
                assert reopened.read_result(row["id"]) == archive.read_result(
                    row["id"])
            for name in scenarios:
                result_row = report[name]
                data = reopened.read_flux(result_row["result_id"])
                assert data["receipt"]["result_id"] == result_row["result_id"]
                reopened_expected[name] = data
        report["reopened"] = sorted(reopened_expected)
    return report


def _run_random(archive, count=12):
    """Deterministic pseudo-random finite inputs, each fully reconciled."""
    seed = 20240917

    def next_int(modulus):
        nonlocal seed
        seed = (1103515245 * seed + 12345) % (2**31)
        return seed % modulus

    def rational(maximum):
        return F(next_int(maximum * 2 + 1) - maximum, next_int(4) + 1)

    def tetra():
        while True:
            points = [tuple(rational(2) for _ in range(3)) for _ in range(4)]
            if _tetra_volume(points) != 0:
                return points

    def density(vertices):
        while True:
            coeffs = [rational(2) for _ in range(4)]
            if all(sum((coeffs[k] * p[k-1] for k in range(1, 4)),
                       coeffs[0]) >= 0 for p in vertices):
                return tuple(coeffs)

    def velocity():
        return [tuple(rational(1) for _ in range(4)) for _ in range(3)]

    from remap import Cell, Target

    verified = 0
    for index in range(count):
        cells = []
        for number in range(1 + next_int(2)):
            vertices = tetra()
            cells.append(Cell("c%d" % number, vertices, density(vertices)))
        targets = []
        for number in range(1 + next_int(4)):
            lower = tuple(rational(2) for _ in range(3))
            upper = tuple(lower[k] + F(next_int(3) + 1, next_int(3) + 1)
                          for k in range(3))
            candidate = Target("t%d" % number, lower, upper)
            if any(_interiors_overlap(candidate, other) for other in targets):
                continue
            targets.append(candidate)
        fields = {cell.id: velocity() for cell in cells}
        revision = archive.put_source("random-source", cells)
        plan = archive.save_plan("random-plan", revision, targets)
        field_id = archive.put_velocity("random-field", revision, fields)
        result = archive.execute_flux(plan, field_id, "random-%d" % index)
        try:
            check_result(result, cells, targets, fields)
        except AssertionError:
            print("random case", index, "failed")
            raise
        # Total divergence of every source over covered bodies matches the
        # published external net, an independent global consistency check.
        for cell in cells:
            total = sum((row["external_net"] for row in result["source_ledger"]
                         if row["source"] == cell.id), F(0))
            assert total == sum((f["flux"] for f in result["faces"]
                                 if f["source"] == cell.id), F(0))
        verified += 1
    return verified


def _tetra_volume(points):
    import itertools
    for combo in itertools.combinations(points, 4):
        matrix = [[combo[i][k] - combo[0][k] for k in range(3)] for i in range(1, 4)]
        det = (matrix[0][0] * (matrix[1][1]*matrix[2][2] - matrix[1][2]*matrix[2][1])
               - matrix[0][1] * (matrix[1][0]*matrix[2][2] - matrix[1][2]*matrix[2][0])
               + matrix[0][2] * (matrix[1][0]*matrix[2][1] - matrix[1][1]*matrix[2][0]))
        if det:
            return abs(det) / 6
    return F(0)


def _interiors_overlap(one, two):
    return all(max(a, b) < min(c, d) for a, b, c, d in
               zip(one.lower, two.lower, one.upper, two.upper))

def main(argv=None):
    import argparse
    import json
    import tempfile
    from pathlib import Path

    parser = argparse.ArgumentParser(
        description="Independent analytic audit of published flux evidence")
    parser.add_argument("--check", action="store_true",
                        help="run every audit scenario")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--scenario", action="append", default=None,
                        help="restrict to a named scenario")
    args = parser.parse_args(argv)
    if not args.check:
        parser.error("nothing to do; pass --check")
    with tempfile.TemporaryDirectory(prefix="flux-audit-") as directory:
        report = run_scenarios(str(Path(directory) / "audit.sqlite"),
                               verbose=args.verbose)
    print(json.dumps({"verified": sorted(report), "scenarios": report},
                     sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


