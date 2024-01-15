### Simulation Loop
*TODO SHOW WHOLE LOOP WITH EXAMPLE!*

One simulation loop should always consist of the following phases:

1. Updating the Container:
   ```cpp
   auto invalidParticles = autoPas.updateContainer();
   ```
   This call will trigger an update of the container inside AutoPas.
   The returned vector `invalidParticles` consists of the particles that were previously owned by this AutoPas container but have left the boundary of this container, i.e., their current position resides outside the container.

2. Handling the leaving particles
    * Apply boundary conditions on them
    * Potentially send them to other mpi-processes, skip this if MPI is not needed
    * Add them to the containers using
       ```cpp
       autoPas.addParticle(particle)
       ```

3. Handle halo particles:
    * Identify the halo particles by use of AutoPas' iterators and send them in a similar way as the leaving particles.
    * Add the particles as haloParticles using
       ```cpp
       autoPas.addHaloParticle(haloParticle)
       ```

4. Perform an iteratePairwise step.
   ```cpp
   autoPas.iteratePairwise(functor);
   ```

## Related Files and Folders
- AutoPasDecl.h
- Simulation.cpp
