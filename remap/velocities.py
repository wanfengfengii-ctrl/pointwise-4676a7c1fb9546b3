"""Immutable affine velocity fields bound to one source revision."""
import json

from .errors import ValidationError
from .models import database_id, identifier, rational, velocity_fields
from .sources import read_source
from .storage import dumps, fetch_one, transaction


def put_velocity(connection, name, source_revision, fields):
    identifier(name)
    database_id(source_revision)
    normalized = velocity_fields(fields)
    source = read_source(connection, source_revision)
    ids = {cell.id for cell in source["cells"]}
    if set(normalized) != ids:
        raise ValidationError(
            "velocity fields must cover exactly the source cell IDs")
    payload = dumps({"version": 1, "source_revision": source_revision,
                     "fields": {cell_id: [[str(x) for x in component]
                                          for component in components]
                                for cell_id, components in normalized.items()}})
    with transaction(connection):
        return connection.execute(
            "INSERT INTO velocity_revisions(name,source_id,payload) VALUES(?,?,?)",
            (name, source_revision, payload)).lastrowid


def read_velocity(connection, revision):
    name, source, payload = fetch_one(
        connection, "SELECT name,source_id,payload FROM velocity_revisions "
        "WHERE id=?", (revision,), "velocity revision")
    data = json.loads(payload)
    fields = {cell_id: tuple(tuple(rational(x) for x in component)
                             for component in components)
              for cell_id, components in data["fields"].items()}
    return {"id": revision, "name": name, "source_revision": source,
            "fields": fields}
