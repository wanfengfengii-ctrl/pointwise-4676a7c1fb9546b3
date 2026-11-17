"""Public flux execution flow: replay, validation, settlement and publish."""
from .errors import ValidationError
from .flux import settle
from .models import database_id, identifier
from .plans import read_plan
from .sources import read_source
from .velocities import read_velocity


def calculate(cells, targets, velocity):
    density = {cell.id: cell.density for cell in cells}
    return settle(cells, targets, density, velocity["fields"])


def execute_flux(archive, plan_id, velocity_revision, request_id):
    identifier(request_id)
    database_id(plan_id)
    database_id(velocity_revision)
    replay = archive._flux_replay(request_id, plan_id, velocity_revision)
    if replay is not None:
        return replay
    plan = read_plan(archive._connection, plan_id)
    velocity = read_velocity(archive._connection, velocity_revision)
    source = read_source(archive._connection, plan["source_revision"])
    if velocity["source_revision"] != plan["source_revision"]:
        raise ValidationError("velocity field and plan use different sources")
    if set(velocity["fields"]) != {cell.id for cell in source["cells"]}:
        raise ValidationError("velocity field does not cover the source cells")
    result = calculate(source["cells"], plan["targets"], velocity)
    result.update({"plan_id": plan_id,
                   "source_revision": plan["source_revision"],
                   "target_revision": plan["target_revision"],
                   "velocity_revision": velocity_revision})
    return archive._publish_flux(plan_id, velocity_revision, request_id, result)
