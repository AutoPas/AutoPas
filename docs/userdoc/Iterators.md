# Iterating Through Particles
Iterators to iterate over particle are provided.
The particle can be accesses using the dereference operators `iter->` or `*iter`.
When created inside a OpenMP parallel region, work is automatically spread over all threads.
```cpp
AUTOPAS_OPENMP(parallel)
for(auto iter = autoPas.begin(); iter.isValid(); ++iter) {
  // user code:
  auto position = iter->getR();
}
```

For convenience the `end()` method is also implemented for the AutoPas class, so you can also use range-based for loops:
```cpp
AUTOPAS_OPENMP(parallel)
for(auto& particle : autoPas) {
  // user code:
  auto position = particle.getR();
}
```

To iterate over a subset of particles, the `getRegionIterator(lowCorner, highCorner)` method can be used:
```cpp
AUTOPAS_OPENMP(parallel)
for(auto iter = autoPas.getRegionIterator(lowCorner, highCorner); iter != autoPas.end(); ++iter) {
  // user code:
  auto position = iter->getR();
}
```

Both `begin()` and `getRegionIterator()` can also take the additional parameter `IteratorBehavior`, which indicates over which particles the iteration should be performed. See [autopas::IteratorBehavior](https://autopas.github.io/doxygen_documentation/git-master/classautopas_1_1options_1_1IteratorBehavior.html) for possible options and details.
The default parameter is `ownedOrHalo`, which is also used for range-based for loops.

Analogously to `begin()`, `cbegin()` is also defined, which guarantees to return a `const_iterator`.

Iterators are not guaranteed to be valid after particle insertion (see [Issue #766](https://github.com/AutoPas/AutoPas/issues/766) for details).
However, particle deletion while iterating is supported via `autoPas.deleteParticle(iterator)`. 
After deletion the `++` operator has to be called:
```cpp
AUTOPAS_OPENMP(parallel)
for(auto iter = autoPas.getIterator(); iter != autoPas.end(); ++iter) {
  autoPas.deleteParticle(iterator);
}
```