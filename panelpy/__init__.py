"""panelpy - metodo dos paineis 2D baseado exclusivamente em escoamento potencial."""
from . import ahmed, analytic, elementary, geometry
from .panels import Panels, influence
from .solver import Solution, rotate, solve

__all__ = ["ahmed", "analytic", "elementary", "geometry", "Panels", "influence",
           "Solution", "solve", "rotate"]
