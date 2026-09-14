"""Four-body transfer approximation; see docs/interplanetary_model.md.

GMs: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc
Planet mean radii: https://ssd.jpl.nasa.gov/planets/phys_par.html
Sun radius is the nominal 695700 km radius used by the existing project.
Jupiter/Saturn use system GMs with DE432s system barycentres, and approximate
planet radii for coarse safety screening. Earth uses Earth-centre GM.
"""

from .constants import G
from .models import CelestialBody


def transfer_bodies():
    # Store equivalent masses so the existing model recovers the reference GM.
    parameters = [
        ("Sun", 1.3271244004193938e11, 695700.0),
        ("Earth", 3.9860043543609598e5, 6371.0084),
        ("Jupiter", 1.2671276480000021e8, 69911.0),
        ("Saturn", 3.7940585200000003e7, 58232.0),
    ]
    return tuple(CelestialBody(name, gm/G, radius) for name, gm, radius in parameters)
