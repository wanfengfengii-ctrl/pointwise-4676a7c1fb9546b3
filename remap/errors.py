class RemapError(Exception):
    """Base class for rejected business operations."""


class ValidationError(RemapError, ValueError):
    pass


class PlanError(RemapError, ValueError):
    pass


class NotFoundError(RemapError, LookupError):
    pass


class RequestConflict(RemapError):
    pass


class BalanceError(RemapError):
    pass
