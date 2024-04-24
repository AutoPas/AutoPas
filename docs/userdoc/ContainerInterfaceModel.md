# Internal Verlet-like container behavior

The behavior described in this section is opaque to users of AutoPas.

We provide a linked cells-like interface to a user of AutoPas, i.e., a container appears to be updated every time step, leaving particles are returned at every time step and particles can be deleted and added independently to the internal state of the container.
Internally we make this possible, by using partial container updates, which collect leaving particles while marking left particles and halo particles as dummy.
Additionally, we maintain a particle buffer that allows to add particles to AutoPas without modifying the underlying container.
This particle buffer is considered in the force calculation and when iterating through particles.
The neighbour lists are updated when a particle moves more than half of the Verlet skin or at the rebuild frequency.

Another performance optimization is made possible by allowing to reuse the neighbor list entries of halo particles of previous time steps.
While the actual particles have already been implicitly deleted (marked as dummy), they still exist.
For their reuse, we try to add halo particles in their original memory location.
If that is, however, not possible, we add them to another particle buffer (the haloParticleBuffer).

Additional information can be found in [PR 642](https://github.com/AutoPas/AutoPas/pull/642)

## Inserting additional particles
Additional particles (e.g. through a grand-canonical thermostat), can be inserted at any point in the simulation loop.
For periodic boundary conditions, or in an MPI-parallel simulation, you, as the user, are responsible for inserting the appropriate halo particles.

## Related Files and Folders
- LeavingParticleCollector.h
- LogicHandler.h
- ParticleContainerInterface.h
