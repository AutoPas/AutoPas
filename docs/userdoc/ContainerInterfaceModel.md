# Container Interface Model

The behavior described in this section is completely hidden from AutoPas users.
To ensure functionality, particles must not move more than `skin / 2` within the specified Verlet rebuild frequency.
This was initially implemented in [PR 642](https://github.com/AutoPas/AutoPas/pull/642), so more details can be found there.

## External Linked Cells-like interface
For an AutoPas user from the outside, it appears that the particle container is fully updated every time they call `AutoPas::updateContainer`.
Particles leaving the domain are always returned, and particles can be added and deleted at any time.
For periodic boundary conditions or in an MPI-parallel simulation, the user is responsible for inserting the appropriate halo particles.

## Internal Verlet-like container behavior
For Verlet list-based containers to perform efficiently, this is a problem, because they rely on their list references to not change until the next list rebuild.
So our solution is to internally use a Verlet-like behavior,  where the actual container is not updated in every time step, leading to particles not necessarily being sorted into their new cells.
We achieve this by avoiding container data structure changes during all updates that do not involve a potential rebuild of neighbor lists.
This means that particles that would normally be deleted, like those leaving the domain or halos, are only marked for deletion (`OwnershipState::dummy`).
They are only really removed during a container update at the end of a rebuild interval or container change.

### `LogicHandler` Buffers
The `LogicHandler` maintains particle buffers that allow the addition of particles to AutoPas without modifying the underlying container.
These buffers are taken into account for the pairwise iteration as well as any other iterators.
There is one pair of buffers for `owned` and `halo` particles per thread to allow parallel particle insertion.

### Reinsert Halos
As mentioned above, halo particles are not deleted immediately but only marked as `dummy`.
Should the halo particle be reinserted before the container structure is updated, we search the vicinity of its coordinates.
If there is a `dummy` particle with the same `id`, we replace it with the new particle.
If nothing is found, the particle is added to a `LogicHandler` halo buffer.
This effectively updates halo particles instead of deleting and reinserting them, making our data structures more stable.

#### Weird Edge Cases
It is still possible that there is a dummy that is not updated, and a new particle is added. Consider the following chain of events:
- The container is `LinkedCells` (or anything based on `CellBlock3D.h`), and a particle is stored in a cell right before the halo region.
- The particle's position changes to the halo region before a data structure update.
- The particle is deleted (= turned `dummy`) by `AutoPas::updateContainer` because it left the domain.
- The particle is immediately added as a halo again. Since halo particles are not allowed to exist in non-halo cells, the existing dummy is not updated, and a new halo particle is added to buffers.

## Related Files and Folders
- LeavingParticleCollector.h
- LogicHandler.h
- ParticleContainerInterface.h
