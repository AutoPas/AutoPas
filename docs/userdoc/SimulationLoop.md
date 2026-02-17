# Simulation Loop
This page describes the assumptions AutoPas makes about how particle simulations look and how AutoPas fits into them.

## High-Level View
AutoPas assumes that any (short-range) particle simulation boils down to three high-level steps:
  1. Particle movement
  2. Particle interaction
  3. Measurements

## Code Example
A conceptual example of a primary simulation loop:
```c++
// main simulation loop
while (needsMoreIterations()) {
   // 1. Particles movement as a result of interactions or simulation forces.
   updatePositions();

   // 1.1 Data structure update
   // Due to AutoPas-external alteration of the particles, AutoPas-internal data structures have to be updated.
   auto emigrants = autopas->updateContainer(); // deletes halo particles, returns particles that are now outside the container

   // 1.2 Particle exchange
   auto [immigrants, haloParticles] = applyBoundaryConditions(emigrants); // user code. Potentially exchanges particles via MPI
   autopas->addParticles(immigrants);           // particles that go into the container
   autopas->addHaloParticles(haloParticles);    // particles that are so close to the container that they have an influence
                                                // on particles inside (periodic boundaries, MPI rank interfaces, ...)

   // 2. Particle Interactions
   YourFunctor functor();                       // User code
   autopas->computeInteractions(&functor);      // Tuning and parallelization happen here

   // 3. Your science goes here.
}
```

## Steps

### 1. Particle Movement
This step is highly user-dependent.
Based on the simulation at hand, particles move due to the result of interaction forces or forces coming from the simulation model.
Nevertheless, positions are usually updated through regular AutoPas iterators
(see [`Iterators.md`](https://github.com/AutoPas/AutoPas/blob/master/docs/userdoc/Iterators.md) for details):
```c++
void updatePositions() {
   for (auto &p : autopas) {
      p.setR(calculateNewPosition(p);
   }
}
```

#### 1.1 Data Structure Update
```c++
auto emigrants = autopas->updateContainer();
```
This method has to be called in every iteration after particles are moved.
It achieves three things (see [`ContainerInterfaceModel.md`](https://github.com/AutoPas/AutoPas/blob/master/docs/userdoc/ContainerInterfaceModel.md) for more details):
  - Internal data structures are updated if necessary.
  - Halo particles are (marked as) deleted.
  - Particles leaving the bounding box of the container are (marked as) deleted and returned as emigrants.

#### 1.2 Particle Exchange
```c++
autopas->addParticles(immigrants);
autopas->addHaloParticles(haloParticles);
```
Depending on the simulation, the user might want to remove and/or insert particles in the local simulation.
Examples of this could be:
  - Particles crossing periodic boundaries (remove in old subdomain, insert in new one).
  - Particles moving from one MPI rank to another.
  - Insertion of particles due to breakup events.
  - Deletion of particles due to outflow condition.

### 2. Particle Interaction
```c++
YourFunctor functor();
autopas->computeInteractions(&functor);
```
Apply your custom force interaction via one or several functors.
This step triggers the tuning.
It internally applies OpenMP parallelization, so do not use it in other thread-parallel environments.

### 3. Measurements
This part encompasses everything else, like sampling particle data, computing macroscopic properties, statistical analysis, ...
For performance reasons, it might be advisable to not do any pairwise particle measurements here, like, for example, pairwise distances for a radial distribution function.
Such `O(N^2)` operations are better handled within or as a separate functor in step 2.
