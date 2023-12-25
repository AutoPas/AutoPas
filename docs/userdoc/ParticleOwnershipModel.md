# Particle Ownership Model
Particles saved in an AutoPas container can be in one of the following states:
* **owned**: Particles that belong to this AutoPas instance.
  These particles are typically inside the boundary of the AutoPas instance.
  If a particle is added via `addParticle()`, it is automatically added as an owned particle.
  An owned particle can explicitly be removed by deleting the particle using an iterator (`autoPas.deleteParticle(iterator)`).
  On an update of the AutoPas container (using `updateContainer()`) owned particles that moved outside the boundary of its parent AutoPas container are returned.
* **halo**: Particles that do not belong to the current AutoPas instance.
  These normally are ghost particles arising from either periodic boundary conditions or particles of a neighboring AutoPas object
  (if you split the entire domain over multiple AutoPas objects, i.e., you use a domain decomposition algorithm).
  The halo particles are needed for the correct calculation of the pairwise forces.
  On update of the AutoPas container, halo particles are deleted (see [Simulation Loop](TODO).
* **dummy**: Particles that are deleted or that act as filler for certain algorithms. These particles do not affect the force calculation.
