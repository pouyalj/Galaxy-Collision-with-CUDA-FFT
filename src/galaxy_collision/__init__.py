"""Galaxy Collision — cross-platform Particle-Mesh N-body simulation.

A modernization of the 2020 PHYD57 project "Galaxy Collisions With CUDA and FFT".
See AGENT.md for the architecture and staged implementation plan.

Stage 0 added the config layer and a trivial ``hello-sim`` smoke run. Stage 1 adds
the unit system (``units``) and the SoA particle/grid data model (``data``). The
physics (IC, deposit, solver, integrator) lands in later stages.
"""

from galaxy_collision import data, units
from galaxy_collision.config import SimConfig, load_config
from galaxy_collision.data import GridState, MemoryEstimate, ParticleState, estimate_memory
from galaxy_collision.units import G

__version__ = "0.0.0"

__all__ = [
    "SimConfig",
    "load_config",
    "units",
    "data",
    "G",
    "estimate_memory",
    "MemoryEstimate",
    "ParticleState",
    "GridState",
    "__version__",
]
