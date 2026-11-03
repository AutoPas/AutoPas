"""Python interface of AutoPas, a node-level auto-tuned particle simulation library.

The central class is :class:`AutoPas`. Its configuration is either a dictionary or a YAML file of
md-flexible, and calling the object starts the simulation or resumes it from its current state::

    import autopas as ap

    sim = ap.AutoPas(
        {
            "cutoff": 2.5,
            "deltaT": 0.002,
            "iterations": 100,
            "Sites": {0: {"mass": 1.0, "epsilon": 1.0, "sigma": 1.0}},
            "Objects": [ap.generators.CubeGrid([10, 10, 10], particle_spacing=1.2)],
        },
        functors=ap.LennardJonesHWY(),
    )
    sim()      # the 100 iterations of the configuration
    sim(50)    # 50 more

See :mod:`autopas.plotting` for plots of the state of a simulation. New functors are written by
deriving from :class:`PairwiseFunctor` or :class:`TriwiseFunctor`.
"""

from ._core import *  # noqa: F403
from ._core import __calculates_globals__, __version__, generators  # noqa: F401
