# Particle Ownership Model

In the spirit of the [container interface model](https://github.com/AutoPas/AutoPas/blob/master/docs/userdoc/ContainerInterfaceModel.md) AutoPas needs a way to distinguish owned and halo particles.
It uses a flagging mechanism, because, depending on the currently used container, they might be mixed in the same data structure.

## Types

Particles in AutoPas have the property `ParticleBase::_ownershipState`, which indicates one of three possible states:

* **owned**: Particles that belong to this AutoPas instance.

  These are the particles one intuitively expects in the simulation.
  They are located inside the boundary of the AutoPas instance.
  On any update of the AutoPas container using `updateContainer()`, owned particles that moved outside the boundary of its parent AutoPas container are turned to `dummy` and returned.

* **halo**: Particles not belonging to the current AutoPas instance.

  These halo particles are needed for the correct calculation of the pairwise forces.
  They usually are ghost particles from periodic boundary conditions, or MPI neighbor ranks.
  On update of the AutoPas container, all halo particles are turned to dummy or deleted (see [Simulation Loop](https://github.com/AutoPas/AutoPas/blob/master/docs/userdoc/SimulationLoop.md)), depending on the underlying data structure.

* **dummy**: Particles that are deleted or that act as filler for specific algorithms.

  These particles do not affect the interaction calculation and are removed from the data container when `updateContainer()` is called, and the data structure is really updated.

## Enforcement

- The ownership state is enforced according to the called function when inserting particles via `addParticle()` or `addHaloParticle()`.
- Due to its close relation to `IteratorBehavior`, static checks are in place to ensure that the binary values for `owned` and `halo` match.
  As explained in the respective classes, this is not easily possible for `dummy`.
- Never change the ownership state of a particle if you do not know precisely what you are doing.
  Especially deletion should always be done through methods of the main interface to avoid confusing internal counters.
