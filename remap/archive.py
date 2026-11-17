import json
from fractions import Fraction
from . import plans, sources, velocities
from .errors import RequestConflict
from .models import database_id
from .storage import connect, dumps, fetch_one, transaction


def _encode(value):
    if isinstance(value, Fraction):
        return {"$r": str(value)}
    if isinstance(value, dict):
        return {k: _encode(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_encode(v) for v in value]
    return value


def _decode(value):
    if isinstance(value, dict):
        if set(value) == {"$r"}:
            return Fraction(value["$r"])
        return {k: _decode(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_decode(v) for v in value]
    return value


class Archive:
    """One local SQLite connection; use a separate instance per worker."""

    def __init__(self, path):
        self._connection = connect(path)

    def close(self):
        self._connection.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def put_source(self, name, cells):
        return sources.put_source(self._connection, name, cells)

    def source_head(self, name):
        return sources.source_head(self._connection, name)

    def read_source(self, revision):
        return sources.read_source(self._connection, revision)

    def list_sources(self):
        return sources.list_sources(self._connection)

    def save_plan(self, name, source_revision, targets):
        return plans.save_plan(self._connection, name, source_revision, targets)

    def read_plan(self, plan_id):
        return plans.read_plan(self._connection, plan_id)

    def list_plans(self):
        return plans.list_plans(self._connection)

    def execute_plan(self, plan_id, request_id, *, observer=None):
        from .settlement import execute_plan
        return execute_plan(self, plan_id, request_id, observer)

    def put_velocity(self, name, source_revision, fields):
        return velocities.put_velocity(self._connection, name, source_revision,
                                       fields)

    def execute_flux(self, plan_id, velocity_revision, request_id):
        from .flux_settlement import execute_flux
        return execute_flux(self, plan_id, velocity_revision, request_id)

    def read_flux(self, result_id):
        database_id(result_id)
        velocity, plan, payload, request = fetch_one(
            self._connection,
            "SELECT f.velocity_id,f.plan_id,f.payload,c.request_id "
            "FROM flux_results f JOIN flux_receipts c ON c.result_id=f.id "
            "WHERE f.id=?", (result_id,), "flux result")
        result = _decode(json.loads(payload))
        result["id"] = result_id
        result["receipt"] = {"request_id": request, "plan_id": plan,
                             "velocity_revision": velocity,
                             "result_id": result_id}
        return result

    def read_result(self, result_id):
        plan, payload, request = fetch_one(
            self._connection,
            "SELECT r.plan_id,r.payload,c.request_id FROM results r "
            "JOIN receipts c ON c.result_id=r.id WHERE r.id=?",
            (result_id,), "result")
        result = _decode(json.loads(payload))
        result["id"] = result_id
        result["receipt"] = {"request_id": request, "plan_id": plan,
                             "result_id": result_id}
        return result

    def list_results(self):
        return [{"id": i, "plan_id": p, "request_id": q} for i, p, q in
                self._connection.execute(
                    "SELECT r.id,r.plan_id,c.request_id FROM results r "
                    "JOIN receipts c ON r.id=c.result_id ORDER BY r.id")]

    def _replay(self, request_id, plan_id):
        row = self._connection.execute(
            "SELECT plan_id,result_id FROM receipts WHERE request_id=?",
            (request_id,)).fetchone()
        if row is None:
            return None
        if row[0] != plan_id:
            raise RequestConflict("request ID already bound to another plan")
        return self.read_result(row[1])

    def _flux_replay(self, request_id, plan_id, velocity_revision):
        row = self._connection.execute(
            "SELECT plan_id,velocity_id,result_id FROM flux_receipts "
            "WHERE request_id=?", (request_id,)).fetchone()
        if row is None:
            return None
        if row[0] != plan_id or row[1] != velocity_revision:
            raise RequestConflict(
                "request ID already bound to another flux execution")
        return self.read_flux(row[2])

    def _publish_flux(self, plan_id, velocity_revision, request_id, result):
        with transaction(self._connection):
            replay = self._flux_replay(request_id, plan_id, velocity_revision)
            if replay is not None:
                return replay
            result_id = self._connection.execute(
                "INSERT INTO flux_results(plan_id,velocity_id,payload) "
                "VALUES(?,?,?)", (plan_id, velocity_revision,
                                  dumps(_encode(result)))).lastrowid
            self._connection.execute(
                "INSERT INTO flux_receipts(request_id,plan_id,velocity_id,"
                "result_id) VALUES(?,?,?,?)",
                (request_id, plan_id, velocity_revision, result_id))
        return self.read_flux(result_id)

    def _publish(self, plan_id, request_id, result):
        with transaction(self._connection):
            replay = self._replay(request_id, plan_id)
            if replay is not None:
                return replay
            result_id = self._connection.execute(
                "INSERT INTO results(plan_id,payload) VALUES(?,?)",
                (plan_id, dumps(_encode(result)))).lastrowid
            self._connection.execute(
                "INSERT INTO receipts(request_id,plan_id,result_id) VALUES(?,?,?)",
                (request_id, plan_id, result_id))
        return self.read_result(result_id)
