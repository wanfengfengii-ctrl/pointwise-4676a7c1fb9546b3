from fractions import Fraction
from .errors import BalanceError
from .geometry import intersection_moments, tetrahedron_moments
from .models import identifier
from .plans import read_plan
from .sources import read_source


def calculate(cells, targets, observer=None):
    transfers, ledger = [], []
    totals = {t.id: Fraction(0) for t in targets}
    for cell in cells:
        full = tetrahedron_moments(cell.vertices)
        original = full.integrate(cell.density)
        covered_mass, covered_volume = Fraction(0), Fraction(0)
        for target in targets:
            overlap = intersection_moments(cell.vertices, target.lower, target.upper)
            if observer is not None:
                observer({"event": "overlap", "source": cell.id, "target": target.id,
                          "moments": overlap})
            mass = overlap.integrate(cell.density)
            if overlap.volume < 0 or overlap.volume > full.volume or mass < 0:
                raise BalanceError("invalid overlap for " + cell.id)
            covered_mass += mass
            covered_volume += overlap.volume
            totals[target.id] += mass
            transfers.append({"source": cell.id, "target": target.id,
                              "volume": overlap.volume, "first": list(overlap.first),
                              "mass": mass})
        remainder = original - covered_mass
        if remainder < 0 or covered_volume > full.volume:
            raise BalanceError("coverage exceeds source " + cell.id)
        ledger.append({"source": cell.id, "total_mass": original,
                       "covered_mass": covered_mass, "uncovered_mass": remainder})
        if sum((r["mass"] for r in transfers if r["source"] == cell.id),
               Fraction(0)) + remainder != original:
            raise BalanceError("source mass does not balance")
    return {"transfers": transfers, "source_ledger": ledger,
            "targets": [{"target": t.id, "mass": totals[t.id]} for t in targets]}


def execute_plan(archive, plan_id, request_id, observer=None):
    identifier(request_id)
    replay = archive._replay(request_id, plan_id)
    if replay is not None:
        return replay
    plan = read_plan(archive._connection, plan_id)
    source = read_source(archive._connection, plan["source_revision"])
    result = calculate(source["cells"], plan["targets"], observer)
    result.update({"plan_id": plan_id, "source_revision": plan["source_revision"],
                   "target_revision": plan["target_revision"]})
    return archive._publish(plan_id, request_id, result)
