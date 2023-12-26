# Custom Applications
AutoPas can be used for arbitrary particle simulations.
The repository comes with a few example applications that can be found in [`applicationLibrary`](https://github.com/AutoPas/AutoPas/blob/master/applicationLibrary). 
To create a new custom application, custom particle and functor classes have to be created. 

## Custom Particles
First you will need to define a particle class which will be passed to AutoPas as template Argument.
For that we provide some basic Particle classes defined in [`molecularDynamicsLibrary`](https://github.com/AutoPas/AutoPas/blob/master/applicationLibrary/molecularDynamics/molecularDynamicsLibrary) or [`SPHLibrary`](https://github.com/AutoPas/AutoPas/blob/master/applicationLibrary/sph/SPHLibrary) that you can use either directly or you can write your own Particle class by inheriting from one of the provided classes or from [`autopas::ParticleBase`](https://github.com/AutoPas/AutoPas/blob/master/src/autopas/particles/ParticleBase.h).

Important parts to implement:
* `enum AttributeNames`
* Definition of a matching `SoAArraysType`
* Getter and setter connecting the `AttributeNames` and actual members.

## Custom Functors
Once you have defined your particle you can start with the functor class.

### Definition
Important parts to implement:
* Actual force calculations: `AoSFunctor()` and all Versions of `SoAFunctor*()`
* Newton3 characteristics of the force: `allowsNewton3()`, `allowsNonNewton3()`
* The calculation of the globals (potential energy, virial) must be implemented in a way so functor calls with newton3 enabled and newton3 disabled within one iteration are possible.
* Input and output variables of the force calculation via: `getComputedAttr()` and `getNeededAttr()`

### Usage
Each functor is applied to AutoPas via:
```cpp
autoPas.iteratePairwise(&myFunctor);
```

### Using multiple functors
AutoPas is able to work with simulation setups using multiple functors that describe different forces.
A good demonstration for that is the sph example found under examples/sph or examples/sph-mpi.
There exist some things you have to be careful about when using multiple functors:
* If you use multiple functors it is necessary that all functors support the same newton3 options.
  If there is one functor not supporting newton3, you have to disable newton3 support for AutoPas by calling
  ```cpp
  autoPas.setAllowedNewton3Options({false});
  ```
* If you have `n` functors within one iteration and update the particle position only at the end or start of the iteration,
  the rebuildFrequency and the samplingRate have to be a multiple of `n`.
* Functors must be marked as (not) relevant for tuning by specifying `Functor::isRelevantForTuning()`.
  Functors marked as relevant should have a near-identical performance profile otherwise the sampling of configurations will be distorted.
  It is recommended, to only mark the most expensive functor as relevant.

## Related Files and Folders
- Functor.h
- ParticleBase.h