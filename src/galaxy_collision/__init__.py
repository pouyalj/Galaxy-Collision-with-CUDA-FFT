"""Galaxy Collision — cross-platform Particle-Mesh N-body simulation.

A modernization of the 2020 PHYD57 project "Galaxy Collisions With CUDA and FFT".
See AGENT.md for the architecture and staged implementation plan.

Stage 0 added the config layer and a trivial ``hello-sim`` smoke run. Stage 1 added the
unit system (``units``) and the SoA particle/grid data model (``data``). Stage 2 added
two-galaxy initial conditions (``ic``). Stage 3 adds the physics: CIC deposit + grad +
gather (``deposit``), the pluggable Poisson solvers (``solver``), the KDK leapfrog
(``integrator``; PM forces are grid-softened, with explicit Plummer softening on the
direct path), and conservation diagnostics (``diagnostics``). The real step loop is wired
into ``sim`` (the ``galaxy-sim`` CLI) and validated by the exit-gate equilibrium tests.
"""

from galaxy_collision import data, ic, units
from galaxy_collision.config import SimConfig, load_config
from galaxy_collision.data import GridState, MemoryEstimate, ParticleState, estimate_memory
from galaxy_collision.ic import GalaxyModel, ICResult, build_ic
from galaxy_collision.units import G

__version__ = "0.0.0"

__all__ = [
    "SimConfig",
    "load_config",
    "units",
    "data",
    "ic",
    "G",
    "estimate_memory",
    "MemoryEstimate",
    "ParticleState",
    "GridState",
    "GalaxyModel",
    "ICResult",
    "build_ic",
    "__version__",
]
