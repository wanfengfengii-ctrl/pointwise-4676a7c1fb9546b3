from fractions import Fraction as F
from remap import Cell, Target

UNIT = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))


def cell(id="s", density=(0, 1, 0, 0), vertices=UNIT):
    return Cell(id, vertices, density)


def whole(id="all"):
    return Target(id, (-1, -1, -1), (2, 2, 2))


def half(id="west"):
    return Target(id, (0, 0, 0), (F(1, 2), 1, 1))


def configured(archive, targets=None, density=(0, 1, 0, 0)):
    revision = archive.put_source("observations", [cell(density=density)])
    plan = archive.save_plan("layout", revision, targets or [half()])
    return revision, plan
