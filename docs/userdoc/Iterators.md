# Iterating Through Particles

Since particles could be stored in many different ways within AutoPas, it is generally not feasible to provide random access to individual particles.
To access particles, iterators are used.
AutoPas provides three interfaces to do this:
1. Classic iterators 
   These `ContainerIterator`s are useable in loops where they provide a maximum of flexibility over the loop flow control.
   They support restriction to certain particle types and regions, as well as parallelization and particle deletion but not addition.
2. `forEach`-style functions (Experimental, no parallelization yet!).
   Similar to iterators, but prohibit interfering with flow control.
3. Functor application functions like `AutoPas::computeInteractions()`.
   Applies the given functor to all particle pairs that are within the cutoff distance.
   
This documentation page shall primarily focus on the first and slightly touch upon the second point.
For the third, refer to [`SimulationLoop.md`](https://github.com/AutoPas/AutoPas/blob/master/docs/userdoc/SimulationLoop.md) for more details.

## Container Iterators

### Range-Based For Loops

The easiest way to iterate over all (= owned and halo) particles of an AutoPas instance is with a [range-based for loop](https://en.cppreference.com/w/cpp/language/range-for):

```c++
for(auto& particle : autoPas) {
  // user code:
  auto position = particle.getR();
}
```

### Manual Iterating

If more fine-grained control is needed over the iteration sequence, the iterator can be accessed and advanced manually:

```c++
for(auto iter = autoPas.begin(); iter != autoPas.end(); ++iter) {
  // user code:
  auto position = iter->getR();
}
```
Analogously to `begin()`, `cbegin()` is also defined, which guarantees to return a `const` iterator.

Some further information to be aware of:
- The iterator can only go forward.
  Its design doesn't necessarily prohibit reverse iteration, but it is (currently) not implemented.
- Adding particles during iteration is considered undefined behavior.
  See [Issue #766](https://github.com/AutoPas/AutoPas/issues/766) for details
- There is no guarantee on the order of iteration.
  If the internal data structure is updated or changed, the iteration order can change as well.
- Alternatively to `iter != autoPas.end()`, we can also use `iter.isValid()`.
  This is equivalent, because `AutoPas::end()` returns a `bool` and `ContainerIterator::operator==(bool)` calls `ContainerIterator::isValid()`.

#### Particle Deletion

AutoPas supports the deletion of particles while iterating without invalidating the iterator.
However, this needs to be done through the member function `AutoPas::deleteParticle()`.
This function also leaves the iterator in a state so that `ContainerIterator::operator++()` needs to be called.

```c++
for(auto iter = autoPas.begin(); iter != autoPas.end(); ++iter) {
  autoPas.deleteParticle(iterator);
}
```

### Restricting Particle Types

Iterators can be restricted to only return particles of a given ownership type.
The default is to cover owned and halo particles.
See [ParticleOwnershipModel.md](https://github.com/AutoPas/AutoPas/blob/master/docs/userdoc/ParticleOwnershipModel.md) for particle ownership types.

```c++
// only iterate owned particles
for(auto iter = autoPas.begin(autopas::IteratorBehavior::owned);
    iter != autoPas.end();
    ++iter) {
  // user code
}
```

Depending on the container, this might be more efficient than iterating everything and calling `continue` on undesired ownership types.

### Region Iterators

Iterators can be restricted to only return particles inside of a given region.
Regions are always cuboids defined by their left-lower-frontal and right-upper-back corner.

```c++
// only iterate particles in a given box
const std::array<double, 3> lowCorner = ...;
const std::array<double, 3> highCorner = ...;
for(auto iter = autoPas.getRegionIterator(lowCorner, highCorner);
    iter != autoPas.end();
    ++iter) {
  // user code
}
```

Region Iterators also take an `IteratorBehavior` as an optional third argument.

### Parallelization

Iteration of the AutoPas object can also be done in parallel via OpenMP.
For this, place any of the examples above in a parallel region.

```c++
AUTOPAS_OPENMP(parallel)
for(auto& particle : autoPas) {
  // user code
}
```

Note that this is equal to `#pragma omp parallel` and **NOT** a `parallel for` clause!
AutoPas detects that it is in a parallel region and spawns multiple iterators spread over the container.
All features described above also work in parallel, including particle deletion.

#### Nested Loops

If iterators are nested, inner loops should not spawn multiple iterators to avoid spreading the work too thin.
This can be achieved by adapting the `IteratorBehavior`, which works similar to a bit vector:
```c++
AUTOPAS_OPENMP(parallel)
for(auto& particle : autoPas) {   // spawns N iterators
  for(auto iter = autoPas.begin(autopas::IteratorBehavior::ownedOrHalo | autopas::IteratorBehavior::forceSequential);
      iter != autoPas.end();
      ++iter) {                   // spawns 1 iterator
    // user code
  }
}
```


## `ForEach`

The `forEach`-style methods aim to provide a similar feature set as `ContainerIterator`.
For this, multiple member functions are available:

| Function | Purpose |
|:---------|:--------|
| `forEach(Lambda, IteratorBehavior)` | Sequential iteration over all particles in the container. |
| `forEachParallel(Lambda, IteratorBehavior)` | Parallel iteration over all particles in the container. |
| `forEachInRegion(Lambda, lowCorner, highCorner, IteratorBehavior)` | Sequential iteration over a given box, same as region iterators. |
| `forEachInRegionParallel(Lambda, lowCorner, highCorner, IteratorBehavior)` | Parallel iteration over a given box, same as region iterators. |

**Parallel functions currently do NOT YET implement any parallelization and an equivalent to their sequential counterparts!**

Here, `Lambda` stands for a function of the signature `(Particle &) -> void`.
All of those functions also exist in a `const` version with the otherwise same signature.

```c++
autoPas.forEachInRegionParallel([](auto &particle) {
  // user code (here e.g. calculating newR)
  particle.setR(newR);
}, lowCorner, highCorner, behavior);
```

### Reductions

On top of that, for every `forEach`-style function, there is also a corresponding `reduce[inRegion][Parallel]()` function.
Here, the lambda function has the form `(Particle &, ResultT &) -> void`
```c++
size_t result = 0;
autoPas.reduceInRegionParallel([](auto &particle, auto &accumulator) {
  // user code
  accumulator += particle.getID();
}, result, lowCorner, highCorner, behavior);
```
