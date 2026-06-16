"""Galaxy Collision — cross-platform Particle-Mesh N-body simulation.

A modernization of the 2020 PHYD57 project "Galaxy Collisions With CUDA and FFT".
See AGENT.md for the architecture and staged implementation plan.

This is Stage 0: scaffold & guardrails. Only the configuration layer and a trivial
``hello-sim`` smoke run exist yet; the physics modules land in later stages.
"""

from galaxy_collision.config import SimConfig, load_config

__version__ = "0.0.0"

__all__ = ["SimConfig", "load_config", "__version__"]
