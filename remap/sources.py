import json
from .errors import ValidationError
from .geometry import tetrahedron_moments
from .models import Cell, identifier, normalize_cell, rational
from .storage import dumps, fetch_one, transaction


def put_source(connection, name, cells):
    identifier(name)
    try:
        cells = tuple(normalize_cell(c) for c in cells)
    except TypeError as exc:
        raise ValidationError("expected source cells") from exc
    if not cells or len({c.id for c in cells}) != len(cells):
        raise ValidationError("source must be nonempty with unique cell IDs")
    for cell in cells:
        tetrahedron_moments(cell.vertices)
        a, b, c, d = cell.density
        if any(a+b*x+c*y+d*z < 0 for x, y, z in cell.vertices):
            raise ValidationError("density must be nonnegative at every vertex")
    cells = sorted(cells, key=lambda c: c.id)
    payload = dumps({"version": 1, "cells": [
        {"id": c.id, "vertices": [[str(x) for x in p] for p in c.vertices],
         "density": [str(x) for x in c.density]} for c in cells]})
    with transaction(connection):
        revision = connection.execute(
            "INSERT INTO source_revisions(name,payload) VALUES(?,?)",
            (name, payload)).lastrowid
        connection.execute(
            "INSERT INTO source_heads(name,source_id) VALUES(?,?) "
            "ON CONFLICT(name) DO UPDATE SET source_id=excluded.source_id",
            (name, revision))
    return revision


def read_source(connection, revision):
    name, payload = fetch_one(
        connection, "SELECT name,payload FROM source_revisions WHERE id=?",
        (revision,), "source revision")
    data = json.loads(payload)
    cells = tuple(Cell(c["id"], tuple(tuple(rational(x) for x in p)
                                    for p in c["vertices"]),
                       tuple(rational(x) for x in c["density"]))
                  for c in data["cells"])
    return {"id": revision, "name": name, "cells": cells}


def source_head(connection, name):
    return fetch_one(connection, "SELECT source_id FROM source_heads WHERE name=?",
                     (name,), "source head")[0]


def list_sources(connection):
    return [{"id": i, "name": name, "is_head": bool(head)}
            for i, name, head in connection.execute(
                "SELECT s.id,s.name,s.id=h.source_id FROM source_revisions s "
                "JOIN source_heads h ON s.name=h.name ORDER BY s.id")]
