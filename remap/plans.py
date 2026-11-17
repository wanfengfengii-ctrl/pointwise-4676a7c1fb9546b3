import json
from itertools import combinations
from .errors import NotFoundError, PlanError, ValidationError
from .models import Target, identifier, normalize_target, rational
from .sources import read_source
from .storage import dumps, fetch_one, transaction


def save_plan(connection, name, source_revision, targets):
    try:
        identifier(name)
        read_source(connection, source_revision)
        targets = tuple(normalize_target(t) for t in targets)
        if not targets or len({t.id for t in targets}) != len(targets):
            raise PlanError("targets must be nonempty with unique IDs")
        for left, right in combinations(targets, 2):
            if all(max(a, b) < min(c, d) for a, b, c, d in
                   zip(left.lower, right.lower, left.upper, right.upper)):
                raise PlanError("target interiors overlap")
    except (ValidationError, NotFoundError, TypeError) as exc:
        raise PlanError(str(exc)) from exc
    targets = sorted(targets, key=lambda t: t.id)
    payload = dumps({"version": 1, "targets": [
        {"id": t.id, "lower": [str(x) for x in t.lower],
         "upper": [str(x) for x in t.upper]} for t in targets]})
    with transaction(connection):
        target_revision = connection.execute(
            "INSERT INTO target_revisions(payload) VALUES(?)", (payload,)).lastrowid
        plan = connection.execute(
            "INSERT INTO plans(name,source_id,target_id) VALUES(?,?,?)",
            (name, source_revision, target_revision)).lastrowid
    return plan


def read_plan(connection, plan):
    name, source, target, payload = fetch_one(
        connection, "SELECT p.name,p.source_id,p.target_id,t.payload "
        "FROM plans p JOIN target_revisions t ON p.target_id=t.id WHERE p.id=?",
        (plan,), "plan")
    targets = tuple(Target(t["id"], tuple(rational(x) for x in t["lower"]),
                           tuple(rational(x) for x in t["upper"]))
                    for t in json.loads(payload)["targets"])
    return {"id": plan, "name": name, "source_revision": source,
            "target_revision": target, "targets": targets}


def list_plans(connection):
    return [{"id": i, "name": n, "source_revision": s, "target_revision": t}
            for i, n, s, t in connection.execute(
                "SELECT id,name,source_id,target_id FROM plans ORDER BY id")]
