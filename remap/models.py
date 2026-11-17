from dataclasses import dataclass
from fractions import Fraction
from .errors import ValidationError


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, str, Fraction)):
        raise ValidationError("expected an integer, rational string, or Fraction")
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValidationError("invalid rational") from exc


def identifier(value):
    if not isinstance(value, str) or not value.strip():
        raise ValidationError("expected a nonempty string ID")
    return value


def vector(values):
    try:
        values = tuple(rational(x) for x in values)
    except TypeError as exc:
        raise ValidationError("expected three coordinates") from exc
    if len(values) != 3:
        raise ValidationError("expected three coordinates")
    return values


@dataclass(frozen=True)
class Cell:
    id: str
    vertices: tuple
    density: tuple


@dataclass(frozen=True)
class Target:
    id: str
    lower: tuple
    upper: tuple


@dataclass(frozen=True)
class Moments:
    volume: Fraction
    first: tuple

    def integrate(self, density):
        a, b, c, d = density
        return a * self.volume + sum(
            (q * m for q, m in zip((b, c, d), self.first)), Fraction(0)
        )


ZERO = Moments(Fraction(0), (Fraction(0),) * 3)


def normalize_cell(cell):
    if not isinstance(cell, Cell):
        raise ValidationError("source members must be Cell instances")
    identifier(cell.id)
    try:
        vertices = tuple(vector(p) for p in cell.vertices)
        density = tuple(rational(x) for x in cell.density)
    except TypeError as exc:
        raise ValidationError("invalid cell sequence") from exc
    if len(vertices) != 4 or len(density) != 4:
        raise ValidationError("a cell needs four vertices and four density coefficients")
    return Cell(cell.id, vertices, density)


def normalize_target(target):
    if not isinstance(target, Target):
        raise ValidationError("plan members must be Target instances")
    identifier(target.id)
    lo, hi = vector(target.lower), vector(target.upper)
    if any(a >= b for a, b in zip(lo, hi)):
        raise ValidationError("target bounds must have positive width")
    return Target(target.id, lo, hi)
