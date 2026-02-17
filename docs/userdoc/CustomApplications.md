# Custom Applications

AutoPas can be used as the core of arbitrary (short-range) particle simulations.
The repository comes with a few example applications that can be found in [`applicationLibrary`](https://github.com/AutoPas/AutoPas/blob/master/applicationLibrary). 
To create a new custom application, custom particle and functor classes have to be created. 

## Custom Particles
Your particle class is an object-oriented representation of your particle model.
For compatibility, it should inherit from [`ParticleBase`](https://github.com/AutoPas/AutoPas/blob/master/src/autopas/particles/ParticleBase.h).
This provides basic features like ID, 3D position, 3D force, and an [`OwnershipState`](https://github.com/AutoPas/AutoPas/blob/master/docs/userdoc/ParticleOwnershipModel.md), as well as, amongst others, functionalities to automatically convert the particle into an SoA representation.

For this to work, the new particle type has to define a few things to help the code generation:
- `enum AttributeNames`:
  This enum has to have an entry for (at least) all class members that are used inside the functor.
  The minimal set of values is what is defined in [`ParticleBase`](https://github.com/AutoPas/AutoPas/blob/master/src/autopas/particles/ParticleBase.h).
  Hence, it is best practice to copy this definition in the new class and extend it as necessary.
- `SoAArraysType`:
  This type definition is a type alias for `autopas::utils::SoAType<...>`, which expects as templates the respective types for the aforementioned `AttributeNames`.
- `get()` and `set()` methods templated with the previously defined `AttributeNames`.

As an example see [`MoleculeLJ`](https://github.com/AutoPas/AutoPas/blob/master/applicationLibrary/molecularDynamics/molecularDynamicsLibrary/MoleculeLJ.h).

### Multiple Particle Types
AutoPas only supports a single particle class.
Therefore, if you want to handle multiple particle types e.g. asteroids vs satellites then you have to differentiate between them manually.
This could be achieved via the introduction of an indicator member in your particle class e.g. `typeId` or having the class itself contain the fields for all types.

## Custom Functors
The functor is a class that defines the interaction of particles, also sometimes referred to as the force kernel. 
AutoPas supports pairwise interactions as well as interactions between triplets of particles.
For compatibility, your functor must inherit from either [`PairwiseFunctor`](https://github.com/AutoPas/AutoPas/blob/master/src/autopas/baseFunctors/PairwiseFunctor.h) or [`TriwiseFunctor`](https://github.com/AutoPas/AutoPas/blob/master/src/autopas/baseFunctors/TriwiseFunctor.h).
These classes suggest how to calculate and store the interactions, as well as some properties of the calculation.

The critical elements to implement are:
- `AoSFunctor()`:
  This function defines how particles interact if both are stored in the Array-of-Structs format, which is exactly what was previously defined as particle class.
  This always has to be implemented.
- All versions of `SoAFunctor()`:
  This function defines the same interaction as `AoSFunctor`, but optimized for different data structure layouts.
  It is possible to not implement this and avoid its usage by never allowing `DataLayout::soa`.
- `allowsNewton3()` and `allowsNonNewton3()`:
  Indicator functions to tell AutoPas if the functor supports optimizations using Newton's third law of motion.
  Should the functor calculate global forces, e.g. potential energy or virial, this must be implemented in a way that supports a mixture of all supported Newton3 modes within one iteration.
- `getNeededAttr()` and `getComputedAttr()`:
  Indicator functions that return a `std::array` of all `AttributeNames`, which the functor needs to load from the particle to perform the calculation, as well as which fields are written.
- `isRelevantForTuning()`:
  Indicator function to tell the tuning mechanism if iterations using this functor should be considered or not.
- `getNumFLOPs()` and `getHitRate()`:
  These functions return the number of FLOPs per traversal of the container and the hit-rate (the ratio of distance calculations
  that lead to functor interactions e.g. force contributions.) These functions are only used if `AUTOPAS_LOG_FLOPS` is
  set to `ON`. If unimplemented, these functions return 0, making the statistics produced by the FLOP logger useless, but
  otherwise not affecting the simulation.
- `getName()`:
  This function should return the name of the functor as a `std::string`.

As an example for a pairwise functor see [`SPHCalcDensityFunctor`](https://github.com/AutoPas/AutoPas/blob/master/applicationLibrary/sph/SPHLibrary/SPHCalcDensityFunctor.h).
As an example for a triwise functor see [`AxilrodTellerMutoFunctor`](https://github.com/AutoPas/AutoPas/blob/master/applicationLibrary/molecularDynamics/molecularDynamicsLibrary/AxilrodTellerMutoFunctor.h).

### Using multiple functors
AutoPas is able to work with simulation setups that use multiple functors of the same or different interaction types to describe different forces.

#### Multiple functors of different interaction types
AutoPas supports pairwise and triwise functors which can also be combined to calculate e.g. different force components.
AutoPas tunes pairwise and triwise functors separately, which means that their configurations are independent.
An example input file to set up such a simulation can be found in [`3BodyTest.yaml`](https://github.com/AutoPas/AutoPas/blob/master/examples/md-flexible/input/3BodyTest.yaml).

#### Multiple functors of the same interaction type
A demonstration of this is the [sph example](https://github.com/AutoPas/AutoPas/blob/master/examples/sph/).
There exist some caveats that have to be considered when using multiple functors of the same interaction type (e.g. all pairwise):
* All functors need to support the same Newton3 options.
  If there is one functor not supporting Newton3, you have to disable Newton3 support for AutoPas by calling
  ```c++
  autopas.setAllowedNewton3Options({false});
  ```
  Otherwise, the algorithm selection might choose a configuration with Newton3 and fail to apply the functor that does not support it.
* If you have `n` functors within one iteration and update the particle position only at the end or start of the iteration,
  the Verlet rebuild interval, and the number of collected samples have to be a multiple of `n`.
* Functors must be marked as (not) relevant for tuning by specifying `Functor::isRelevantForTuning()`.
  Functors marked as relevant should have a near-identical performance profile, otherwise, the sampling of configurations will be distorted.
  It is recommended, to only mark the most expensive functor as relevant.
