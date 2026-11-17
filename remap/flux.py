"""Exact directed boundary evidence and per-face mass flux.

Every execution rebuilds the whole boundary from the immutable source,
target and velocity versions; nothing computed for one execution is reused
by another.
"""
from fractions import Fraction
from itertools import combinations

from .errors import BalanceError
from .fluxgeom import (MONOMIALS, clean_polygon, decompose_outside, dot,
                       intersect_axis, polygon_area_vector, polygon_moments)
from .geometry import (_face_hull, _meet, _planes, cross, intersection_moments,
                       sub, tetrahedron_moments)

ZERO3 = (Fraction(0), Fraction(0), Fraction(0))


# --------------------------------------------------------------- plane algebra

def _source_plane(vertices, omitted):
    face = [vertices[i] for i in range(4) if i != omitted]
    p, q, r = face
    normal = cross(sub(q, p), sub(r, p))
    if dot(normal, sub(vertices[omitted], p)) > 0:
        normal = tuple(-x for x in normal)
    return normal, dot(normal, p)


def _box_plane(axis, side, lower, upper):
    n = tuple(Fraction(int(i == axis)) for i in range(3))
    if side == "lower":
        return tuple(-x for x in n), -lower[axis]
    return n, upper[axis]


def _same_halfspace(n1, d1, n2, d2):
    scale = None
    for a, b in zip(n1, n2):
        if a or b:
            if not a or not b:
                return False
            candidate = a / b
            if scale is None:
                scale = candidate
            elif scale != candidate:
                return False
    return scale is not None and scale > 0 and d1 == scale * d2


def _orient(poly, normal):
    poly = clean_polygon(poly)
    if len(poly) < 3:
        return poly
    if dot(polygon_area_vector(poly), normal) < 0:
        return list(reversed(poly))
    return poly


def _merged_planes(vertices, lower, upper, target_id):
    """Collapse coincident source/target planes, retaining both origins."""
    raw = []
    for omitted in range(4):
        normal, bound = _source_plane(vertices, omitted)
        face = sorted(list(vertices[i]) for i in range(4) if i != omitted)
        raw.append((normal, bound, face, None))
    for axis in range(3):
        for side in ("lower", "upper"):
            normal, bound = _box_plane(axis, side, lower, upper)
            raw.append((normal, bound, None, [target_id, axis, side]))
    merged = []
    for normal, bound, source_face, target_face in raw:
        for item in merged:
            if _same_halfspace(normal, bound, item["normal"], item["bound"]):
                if source_face is not None:
                    item["source_faces"].append(source_face)
                if target_face is not None:
                    item["target_faces"].append(target_face)
                break
        else:
            merged.append({
                "normal": normal, "bound": bound,
                "source_faces": [source_face] if source_face is not None else [],
                "target_faces": [list(target_face)] if target_face else []})
    return merged


def body_pieces(vertices, lower, upper, target_id):
    """Facets of the positive-volume tetrahedron/box intersection.

    Empty for a non-positive-volume intersection (edge, point or face contact).
    """
    tetrahedron_moments(vertices)
    planes = _planes(vertices, lower, upper)
    points = set()
    for triple in combinations(planes, 3):
        point = _meet(*triple)
        if point is not None and all(dot(n, point) <= d for n, d in planes):
            points.add(point)
    points = sorted(points)
    if len(points) < 4:
        return []
    if intersection_moments(vertices, lower, upper).volume <= 0:
        return []
    pieces = []
    for item in _merged_planes(vertices, lower, upper, target_id):
        face = [p for p in points if dot(item["normal"], p) == item["bound"]]
        if len(face) < 3:
            continue
        hull = _orient(_face_hull(face, item["normal"]), item["normal"])
        if len(hull) < 3 or polygon_area_vector(hull) == ZERO3:
            continue
        target = item["target_faces"][0] if item["target_faces"] else None
        pieces.append({
            "vertices": hull, "normal": item["normal"],
            "source_faces": [list(map(tuple, f)) for f in item["source_faces"]],
            "target_faces": [list(t) for t in item["target_faces"]],
            "axis": target[1] if target else None,
            "side": target[2] if target else None})
    return pieces


# ------------------------------------------------------------------- flux math

def _quadratic_terms(density, component):
    a, b, c, d = density
    e, f, g, h = component
    return {
        "1": a * e,
        "x": a * f + b * e,
        "y": a * g + c * e,
        "z": a * h + d * e,
        "xx": b * f,
        "xy": b * g + c * f,
        "xz": b * h + d * f,
        "yy": c * g,
        "yz": c * h + d * g,
        "zz": d * h}


def face_flux(moments, density, velocity):
    """Exact integral of ``rho * v . n`` over a facet using oriented moments."""
    total = Fraction(0)
    for axis in range(3):
        terms = _quadratic_terms(density, velocity[axis])
        for monomial in MONOMIALS:
            coefficient = terms[monomial]
            if coefficient:
                total += coefficient * moments[monomial][axis]
    return total


# --------------------------------------------------------------- reconciliation

def _record(source, target, vertices, normal, source_faces, target_faces, mate):
    return {"source": source, "target": target,
            "vertices": [tuple(p) for p in vertices],
            "normal": normal,
            "origins": {"source_faces": [list(map(tuple, f))
                                         for f in source_faces],
                        "target_faces": [list(t) for t in target_faces]},
            "mate": mate}


def _other_side(side):
    return "lower" if side == "upper" else "upper"


def build_evidence(cells, targets):
    """Ordered face records for one immutable input snapshot."""
    bodies = {}
    for cell in cells:
        for target in targets:
            bodies[(cell.id, target.id)] = body_pieces(
                cell.vertices, target.lower, target.upper, target.id)
    faces = []
    for cell in cells:
        for target in targets:
            pieces = bodies[(cell.id, target.id)]
            if not pieces:
                continue
            for piece in pieces:
                if piece["axis"] is None:
                    faces.append(_record(cell.id, target.id, piece["vertices"],
                                         piece["normal"], piece["source_faces"],
                                         piece["target_faces"], None))
                    continue
                axis, side = piece["axis"], piece["side"]
                bound = (target.upper if side == "upper"
                         else target.lower)[axis]
                other_side = _other_side(side)
                holes, pairs = [], []
                for other in targets:
                    if other.id == target.id:
                        continue
                    if side == "upper":
                        if other.lower[axis] != target.upper[axis]:
                            continue
                    else:
                        if other.upper[axis] != target.lower[axis]:
                            continue
                    # The neighbour's own facet, not one of this body's.
                    other_piece = None
                    for candidate in bodies[(cell.id, other.id)]:
                        if (candidate["axis"] == axis
                                and candidate["side"] == other_side):
                            other_piece = candidate
                            break
                    if other_piece is None:
                        continue
                    common = _orient(intersect_axis(
                        piece["vertices"], other_piece["vertices"], axis, bound),
                        piece["normal"])
                    if polygon_area_vector(common) == ZERO3:
                        continue
                    holes.append(common)
                    if target.id < other.id:
                        pairs.append((other, other_piece, common))
                leftovers = decompose_outside(piece["vertices"], axis, bound,
                                              holes)
                for leftover in leftovers:
                    leftover = _orient(leftover, piece["normal"])
                    if polygon_area_vector(leftover) == ZERO3:
                        continue
                    faces.append(_record(cell.id, target.id, leftover,
                                         piece["normal"], piece["source_faces"],
                                         piece["target_faces"], None))
                for other, other_piece, common in pairs:
                    poly_a = _orient(common, piece["normal"])
                    poly_b = list(reversed(poly_a))
                    index = len(faces)
                    faces.append(_record(
                        cell.id, target.id, poly_a, piece["normal"],
                        piece["source_faces"], piece["target_faces"],
                        index + 1))
                    faces.append(_record(
                        cell.id, other.id, poly_b,
                        tuple(-x for x in piece["normal"]),
                        other_piece["source_faces"],
                        other_piece["target_faces"], index))
    return faces


def settle(cells, targets, density_by_id, velocity_by_id):
    """Compute and self-check the complete published flux result body."""
    faces = build_evidence(cells, targets)
    volumes = {(c.id, t.id): intersection_moments(
        c.vertices, t.lower, t.upper).volume for c in cells for t in targets}
    for index, face in enumerate(faces):
        face["id"] = index
        face["moments"] = polygon_moments(face["vertices"])
        face["flux"] = face_flux(face["moments"],
                                 density_by_id[face["source"]],
                                 velocity_by_id[face["source"]])
    by_id = {f["id"]: f for f in faces}
    for face in faces:
        if face["mate"] is None:
            continue
        other = by_id.get(face["mate"])
        if other is None or other["mate"] != face["id"]:
            raise BalanceError("face mate is not symmetric")
        if other["source"] != face["source"]:
            raise BalanceError("paired faces from different sources")
        if other["flux"] != -face["flux"]:
            raise BalanceError("paired face fluxes do not cancel")
    targets_out = {t.id: Fraction(0) for t in targets}
    external = {c.id: Fraction(0) for c in cells}
    for face in faces:
        targets_out[face["target"]] += face["flux"]
        if face["mate"] is None:
            external[face["source"]] += face["flux"]
    for cell in cells:
        for target in targets:
            body = [f for f in faces
                    if f["source"] == cell.id and f["target"] == target.id]
            area = [sum((f["moments"]["1"][k] for f in body), Fraction(0))
                    for k in range(3)]
            divergence = sum((f["moments"]["x"][0] + f["moments"]["y"][1]
                              + f["moments"]["z"][2] for f in body), Fraction(0))
            if area != [Fraction(0)] * 3:
                raise BalanceError("boundary is not closed for " + cell.id)
            if divergence != 3 * volumes[(cell.id, target.id)]:
                raise BalanceError("boundary coverage is incomplete for "
                                   + cell.id)
        net = sum((f["flux"] for f in faces if f["source"] == cell.id),
                  Fraction(0))
        if net != external[cell.id]:
            raise BalanceError("source face flux does not reconcile")
    exchanges = []
    for face in faces:
        if face["mate"] is None:
            continue
        other = by_id[face["mate"]]
        if face["target"] < other["target"]:
            exchanges.append({"source": face["source"],
                              "a": face["target"], "b": other["target"],
                              "a_face": face["id"], "b_face": other["id"],
                              "flux": face["flux"]})
    return {
        "faces": [_public(f) for f in faces],
        "exchanges": exchanges,
        "targets": [{"target": t.id, "net_outward": targets_out[t.id]}
                    for t in targets],
        "source_ledger": [{"source": c.id, "external_net": external[c.id]}
                          for c in cells],
    }


def _public(face):
    return {"id": face["id"], "source": face["source"],
            "target": face["target"],
            "vertices": [list(p) for p in face["vertices"]],
            "origins": {
                "source_faces": [[list(p) for p in f]
                                 for f in face["origins"]["source_faces"]],
                "target_faces": [list(t) for t in
                                 face["origins"]["target_faces"]]},
            "moments": {m: list(face["moments"][m]) for m in MONOMIALS},
            "flux": face["flux"],
            "mate": face["mate"]}
