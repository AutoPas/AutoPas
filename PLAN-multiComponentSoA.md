# Multi-component SoAs for AutoPas (rework of `rigidBodySoALoading`)

> **Self-contained hand-off document.** Written 2026-10-20 to be executed in a fresh session (human or Claude Code)
> on a different machine. It assumes no access to the conversation that produced it. The old branch
> `origin/rigidBodySoALoading` is kept only as a reference: `git show origin/rigidBodySoALoading:<path>`.
>
> **Status: provisional.** PR-1 (`refactor/particleConcept`) is being finished separately, outside this plan.
> **Rewrite this plan once PR-1 has landed.** The sections marked *(depends on PR-1)* are the ones expected to change.

## PR sequence

| PR | Branch (suggested) | Content | Review |
|----|--------------------|---------|--------|
| **PR-1** | `refactor/particleConcept` (exists) | `autopas::utils::ParticleType` concept. Done by the user, separately from this plan. | normal |
| **PR0** | `cpp20SoA` | C++20 modernisation of the *existing* SoA handling. No behaviour change. | normal |
| **PR1** | `multiComponentSoA` | Core multi-component SoA (`src/autopas`, tests, docs) | **after PR2 is complete** |
| **PR2** | `multisiteSoASites` (stacked on PR1) | `MultisiteMoleculeLJ` site components, `LJMultisiteFunctor` SoA port, md-flexible | together with PR1 |

- PR0 is independent of PR-1 and could start in parallel, but expect small conflicts in `Functor.h`, because PR-1
  constrains its template parameter. Rebase onto master after PR-1 merges.
- PR1 is not reviewed until PR2 is complete. **PR1's API may still change while working on PR2, but should be
  mostly final before PR2 starts.**

## Context

AutoPas SoAs are flat: one aligned array per particle attribute (`SoAArraysType = SoAType<Ts...>::Type`, a
`std::tuple` of aligned vectors), filled by the generic `Functor::SoALoader`/`SoAExtractor` from the attribute
lists that each functor declares in `getNeededAttr()`/`getComputedAttr()`.

Multi-site molecules (`mdLib::MultisiteMoleculeLJ`) store only centre-of-mass (CoM), quaternion and `typeId`.
Site positions are derived from these plus the `ParticlePropertiesLibrary` (PPL). Because the SoA cannot hold
per-site data, `LJMultisiteFunctor`'s SoA kernels rebuild site-position, site-type and site-force arrays and
reduce forces and torques **inside every `SoAFunctorSingle`/`Pair`/`Verlet` call**. That is O(N) extra work and
memory traffic per cell pair. See `LJMultisiteFunctor.h`, `SoAFunctorSingle`, from roughly line 570 on master.

**Goal:** let a particle declare variable-count *components* (e.g. LJ sites). The SoA then gets extra
*additional partitions* that `SoALoader` fills **once**. Functors read them like any other SoA arrays, and
`SoAExtractor` reduces them back to the particle once, for example site forces into molecule force and torque.

### Layout (unchanged from the old branch)
- The **main partition** is today's SoA, with one array per attribute and length N (number of particles).
- For each **additional partition type** (component type) there is a `std::vector` of partitions. The vector's size
  is the **depth**, which must be at least the largest number of components of that type on any particle in the SoA.
  The partition at depth `d` holds component `d` of every particle, again with arrays of length N.
- Particles with fewer components get "absent" entries. These are value-initialised, and since
  `OwnershipState{} == dummy` (0), a component partition containing an `ownershipState` attribute masks itself.
- This layout fits AutoPas well. Views, appends, offsets and Verlet-list indices remain plain molecule-index
  ranges, and for a fixed depth the inner loops are contiguous and SIMD-friendly. The cost is wasted slots for
  mixtures of molecules with different site counts. That cost is acceptable and documented.

## What the old branch did, and why we are not merging it

The branch has 43 commits, last touched 2024-11-26, and its merge-base is `09bebc93a`. Upstream master is
**1,685 commits ahead**, and **the branch does not compile**. It did the following:
- Renamed `SoA` → `SoAPartition` and `SoAView` → `SoAPartitionView`, and added a new `SoA` = main partition +
  `tuple<vector<SoAPartition<Ti>>...>`.
- Made every particle declare `AttributeTypes` (a `std::tuple` of types) plus
  `SoAPartitionType<AttributeTypes, attrs...>` subsets, which changed every particle class.
- Added `getNeededAdditionalAttr()`/`getComputedAdditionalAttr()` on functors and rewrote
  `SoALoader`/`SoAExtractor` with index_sequence helper chains.
- Included uncommitted WIP: per-site getters on `MultisiteMoleculeLJ` (`getSitePosX(i)`, ...), a static PPL
  reference, and `rotatePositionX/Y/Z` and `crossX/Y/Z` helpers for per-coordinate getters.

Defects to avoid repeating in the new design:
1. **Partition indexing.** Storage is positional (`std::get<attr>` on the partition tuple), so a subset of a single
   `AttributeNames` enum (e.g. `sitePosX = 26`) indexes out of range. `MultisiteMoleculeLJ`'s main partition also
   dropped `id`, which shifted every index after it.
   → **Each component type gets its own enum, starting at 0.**
2. **`SoAView` allocated a `std::vector<SoAPartitionView>` per view.** Views are created in hot loops and passed
   by value.
   → **Views stay {`SoA*`, start, end} and compute component pointers on demand.**
3. **No depth sizing.** Nothing set the depth from the particles. The skip-resize path, where Verlet-list containers
   load cells in parallel at offsets, cannot grow the depth safely.
4. **Per-coordinate setters cannot compute torque.** τ = r × F needs all three force components at once.
   → **Get and set a whole component at a time.** This also makes `rotatePositionX/Y/Z` and `crossX/Y/Z`
   unnecessary, so do not port them.
5. **Many compile errors.** Examples: `soaTmp->clear()` on values, `for (auto soaTmp : ...)` resizing copies,
   `std::get<i, f()>` syntax, `SoA::begin<Functor,bool>` returning a reference to a temporary, and the Extractor
   sizing its loop by the *needed* attributes.
6. **Breaking change** for every particle in core, applicationLibrary, examples, tests and external users.

## Decisions (agreed with the user)

- **Fresh branches from current `origin/master`.** Port the ideas, not the commits.
- **Non-breaking opt-in.** `SoAArraysType = SoAType<...>::Type` stays valid unchanged. Only multi-component
  particles use a new descriptor type. Existing particles, functors and tests must compile **unchanged** in PR0 and
  PR1.
- **Particle-side site logic.** `MultisiteMoleculeLJ` gets a static PPL pointer, which is set once at init.
  The generic `Functor::SoALoader`/`SoAExtractor` stays particle-agnostic.
- **PR1 is exercised by a toy multi-component test particle** in `tests/testAutopas`, with no dependency on
  applicationLibrary. Everything md-specific goes in PR2.

## PR0: C++20 modernisation of the existing SoA handling (no behaviour change)

Scope: `src/autopas/utils/{SoA,SoAView,SoAStorage,SoAType,SortedSoAView}.h`, `src/autopas/baseFunctors/Functor.h`
(loader/extractor), and `LCC01Traversal.h::appendNeeded`.
- Replace index_sequence `*Impl` helper chains with **templated lambdas + fold expressions** (`[&]<size_t... I>(std::index_sequence<I...>) {...}(std::make_index_sequence<N>{})`):
  - `Functor::SoALoaderImpl`/`SoAExtractorImpl`.
  - `SoA::append_impl`/`appendSingleArray` (both overloads).
  - The recursive `SoA::read_impl`/`write_impl` (end-of-recursion overloads included).
  - `SortedSoAView::packNeededImpl`/`zeroComputedImpl`/`scatterBackImpl`.
  - `LCC01Traversal::appendNeeded`'s index_sequence parameter.
- Inside these lambdas use `Functor_T::getNeededAttr()[I]` directly. Capturing a local constexpr array can trip
  GCC's capture rules.
- Resolve the `/* @todo C++20 [[unlikely]] */` in `SoAView`. Drop redundant `typename` where C++20 allows.
- Constrain template parameters with concepts where it clarifies intent, e.g. attribute packs and `ValueArrayType`.
  Coordinate with PR-1's `ParticleType` so the two don't diverge.
- Leave `SoAStorage::apply`'s "expansion statement" todo; that feature is C++26.
- **Out of scope:** particles' `std::enable_if_t` getters (`ParticleBase`, `MoleculeLJ`, `MultisiteMoleculeLJ`,
  `SPHParticle`, test particles). If those are wanted, they belong with PR-1 or a separate PR.
- Verify with the full test suite unchanged, plus a quick md-flexible LJ SoA timing check showing no regression.
  The loader loop should compile to the same code.

## PR1 target design (core, builds on PR0)

### Types: `src/autopas/utils/SoAType.h`
```cpp
// unchanged
template <typename... soatypes> struct SoAType { using Type = std::tuple<std::vector<soatypes, AlignedAllocator<soatypes>>...>; };

// new: descriptor for particles with additional component partitions
template <class MainArraysType, class... AdditionalArraysTypes>
struct MultiComponentSoAType {};

// traits (+ a concept), plain tuples map to {Main = T, Additional = std::tuple<>}
template <class T> struct SoAArraysTraits;                // ::MainArraysType, ::AdditionalArraysTypes (std::tuple<...>)
template <class T> concept MultiComponentSoAArrays = ...; // true only for MultiComponentSoAType<...>
template <class T> using SoAMainArraysType_t = typename SoAArraysTraits<T>::MainArraysType;
```

### `SoA<SoAArraysType>`: `src/autopas/utils/SoA.h`
- Move today's `SoA` body, as modernised in PR0, into a new internal class `SoAPartition<ArraysType>`
  (new file `SoAPartition.h`). Do this first as a pure-refactor commit with no behaviour change.
- `SoA` holds `SoAPartition<Main> _mainPartition;` and
  `[[no_unique_address]] std::tuple<std::vector<SoAPartition<Ai>>...> _additionalPartitions;`.
  The tuple is `std::tuple<>` for simple particles, so there is **zero overhead**. Assert this with a
  `sizeof` check in the tests.
- Every existing method keeps its signature: `push`, `write`, `read`, `readMultiple`, `writeMultiple`,
  `begin<attr>()` (const and non-const), `append<attrs...>`, `constructView`, `size`.
  - Per-attribute methods forward to the main partition.
  - Whole-SoA operations (`resizeArrays`, `clear`, `swap`, `pop_back`, `append(SoA)`, `append(SoAView)`) apply to
    the main partition and to every additional partition at every depth.
- New methods: `begin<type, attr>(size_t depth)` (const and non-const), `getDepth<type>()`,
  `resizeDepth<type>(n)`, and `getDepths()` (returns `std::array`). New partitions are sized to the current length
  and value-initialised, i.e. absent.
- `append` takes the max of the two depths. Partitions deeper than the other SoA's depth get value-initialised rows.
  The selective `append<attrs...>` used by LCC01 (N3-off combined SoA) stays selective on the main partition and
  appends additional partitions in full. Note the selective version as a possible follow-up.
- Carry over the old branch's class doc comment from `git show origin/rigidBodySoALoading:src/autopas/utils/SoA.h`,
  lines 19-60 (TIP4P example and depth explanation), adapted to the new names.

### `SoAView`: `src/autopas/utils/SoAView.h`
- Keep the layout {`SoA*`, start, end}. Add `begin<type, attr>(depth)` (returns `_soa->begin<type,attr>(depth) +
  _startIndex`) and `getDepth<type>()`. There is no allocation.

### `Functor`: `src/autopas/baseFunctors/Functor.h`
- Add defaults `constexpr static std::tuple<> getNeededAdditionalAttr()` and `getComputedAdditionalAttr()`.
  - A functor for multi-component particles returns
    `std::tuple{std::array{Particle_T::SiteAttributeNames::posX, ...}, ...}`, with one array per additional
    partition type.
  - `static_assert` that the tuple size is 0 ("needs nothing") or equals the number of additional partition types.
- Extend PR0's lambda-based `SoALoader`/`SoAExtractor` for additional partitions, only under
  `if constexpr (MultiComponentParticleType<Particle_T>)`. The main-partition path is unchanged.
  - **Depth sizing (non-skip path):** depth = max components in the cell (`offset == 0`), or the max of the current
    depth and the cell's max (`offset > 0`).
  - **Load:** for each type T, each depth `d < soa.getDepth<T>()` and each particle:
    - If `d < p.getNumComponents<T>()`, then `auto c = p.getComponent<T>(d)` and write the needed attributes via
      `std::get<attr>(c)`.
    - Otherwise write `value_type{}` for the needed attributes.
  - **Extract:** for `d < p.getNumComponents<T>()`, call `p.template setComponent<T, computedAttrsOfT>(d, valuesTuple)`.
    The `std::array` is passed as an NTTP via `template <auto attrs>`.
- **Skip-resize path:** add a public helper on `Functor`, suggested name
  `SoAResizeForCells(const Cells &cells, SoA &soa, size_t length)`. It runs `resizeArrays(length)` and, for
  multi-component particles, a serial max-components pass and `resizeDepth`. `SoALoader` with `skipSoAResize == true`
  then asserts that the depth is already sufficient.
  - Replace `_soa.resizeArrays(offsets.back())` with this helper at every `/*skipSoAResize*/ true` call site.
    On 2024 master these were `VerletNeighborListAsBuild.h`, `VLListIterationTraversal.h` and
    `VLCNeighborListInterface.h`. **Re-grep on current master** for more.

### Particle opt-in contract *(depends on PR-1)*
Express the contract as a concept refining PR-1's `autopas::utils::ParticleType` (`src/autopas/utils/ParticleConcept.h`
on `refactor/particleConcept`). Suggested name `MultiComponentParticleType<P>`: it requires `ParticleType<P>`,
`MultiComponentSoAArrays<typename P::SoAArraysType>`, and the members below. `Functor` and the tests then use the
concept rather than ad-hoc detection. Provisional shape:
```cpp
enum AttributeNames : int { ptr, id, posX, ..., ownershipState };              // main, as today
enum AdditionalPartitionTypeNames : size_t { site };                           // index into additional partitions
enum SiteAttributeNames : int { posX, posY, posZ, forceX, forceY, forceZ, ownershipState }; // own enum, from 0
using SoAArraysType = autopas::utils::MultiComponentSoAType<
    autopas::utils::SoAType<Self *, size_t, double, ...>::Type,               // main
    autopas::utils::SoAType<double, double, double, double, double, double, autopas::OwnershipState>::Type>;
template <AdditionalPartitionTypeNames T> size_t getNumComponents() const;
template <AdditionalPartitionTypeNames T> auto getComponent(size_t depth) const;  // tuple of all attrs, enum order
template <AdditionalPartitionTypeNames T, auto attrs> void setComponent(size_t depth, const auto &values);
// main get/set use SoAMainArraysType_t<SoAArraysType> instead of SoAArraysType
```
The semantics of `setComponent` are defined by the particle. For example, site forces are *accumulated* into the
molecule's force and torque.

When PR-1 has landed, check the following:
- Does `ParticleType` still hold for a particle whose `SoAArraysType` is the descriptor? Its SoA requirements
  (`get<ptr>`, `get<posX>`, ...) should not depend on `SoAArraysType`.
- Is `checkFunctorType.h` affected, or any other template that PR-1 constrained?

### Other consumers to update or guard
- `SortedSoAView.h`: it does not permute additional partitions. Add a `static_assert` (or require
  `supportsSoASorting == false`) for multi-component particles.
- `LCC01Traversal.h` `appendNeeded` and the `LCC04SoACellHandler.h` buffers: verify that they handle depth through
  `SoA::append`/`resizeArrays`.
- `docs/userdoc/CustomApplications.md`: document the opt-in.

### PR1 commit sequence
1. Pure refactor: `SoA` → `SoAPartition` + thin `SoA` wrapper. All tests pass with no test changes.
2. `MultiComponentSoAType`, traits/concepts (incl. `MultiComponentParticleType`), and additional partitions in
   `SoA`, plus unit tests.
3. `SoAView` component access.
4. `Functor` component load/extract + depth sizing.
5. Skip-resize helper + container call sites. Guard `SortedSoAView`. Check LCC01/LCC04.
6. Toy particle `tests/testAutopas/testingHelpers/MultiComponentParticle.h`:
   - Two component types with explicitly stored, per-particle variable counts. No PPL.
   - A toy pairwise functor (component-component interaction with a CoM cutoff).
7. Tests (see Verification), then the doc update.

## PR2 (applicationLibrary + md-flexible), stacked on PR1
- `MultisiteMoleculeLJ`:
  - Opt into `MultiComponentSoAType` with an LJ-site partition: site posX/Y/Z, siteTypeId, site forceX/Y/Z and
    site ownershipState.
  - Add a static `const ParticlePropertiesLibrary<double, size_t> *` plus a setter. md-flexible uses
    `ParticlePropertiesLibrary<FloatPrecision, size_t>` (`examples/md-flexible/src/TypeDefinitions.h`), so check
    the types match.
  - `getNumComponents` = `PPL->getNumSites(typeId)`.
  - `getComponent` = one `rotatePosition(q, sitePos[d]) + r`, the site type, zero force and the molecule's
    ownership state.
  - `setComponent<site,{fX,fY,fZ}>` adds F to `_f` and `cross(rotated, F)` to `_torque`.
  - Watch for an include cycle with `ParticlePropertiesLibrary.h`. Move the bodies to the .cpp if needed.
- `ParticlePropertiesLibrary::getSitePositions`/`getSiteTypes`: return `const &` (from the WIP).
- `LJMultisiteFunctor`:
  - Declare the additional attrs.
  - Rewrite `SoAFunctorSingle`/`Pair`/`Verlet` to loop over molecule i → its sites → depth b → j (contiguous,
    masked by the CoM cutoff and the site ownership state), accumulating into the site force arrays.
  - Remove the per-call site-array construction and the torque post-processing.
  - Drop main force/torque from the SoA needed/computed lists if they are no longer touched.
  - Update the FLOP accounting.
- md-flexible: set the static PPL wherever the PPL is finalised; grep for `calculateMixingCoefficients`.
  Test fixtures that build a PPL must set it too.
- Do not port `MultisiteMolLJExactSitePosVar` or `rotatePosition{X,Y,Z}`/`cross{X,Y,Z}` from the WIP.
- Feed back into PR1: any PR1 API change found here is allowed, because PR1 is not yet reviewed. Keep PR1's
  commits coherent, e.g. with fixup commits on the PR1 branch, then rebase PR2.

## Step 0: hand-off (done on the current laptop before switching machines)
1. On `rigidBodySoALoading`:
   - Copy this file to `PLAN-multiComponentSoA.md` at the repo root.
   - Commit it together with the uncommitted WIP (`MultisiteMoleculeLJ.{h,cpp}`, `ParticlePropertiesLibrary.h`, and
     the staged `MultisiteMolLJExactSitePosVar.{h,cpp}`). Mark the commit message "WIP, does not compile".
   - `git push origin rigidBodySoALoading`. This is a fast-forward, because the remote is at `fef7cc5`.
2. On the workstation:
   - Run `git fetch origin` and read the plan with `git show origin/rigidBodySoALoading:PLAN-multiComponentSoA.md`.
   - Finish PR-1 first. Then **rewrite this plan** (a Claude Code session can do it from the plan plus the merged
     PR-1 diff).
   - Then branch PR0 from `origin/master`.

## Verification
- **All PRs:** GCC and Clang, Release and Debug, with OpenMP ON and OFF. The skip-resize parallel path only matters
  with OpenMP. Run the formatting targets (`AUTOPAS_FORMATTING_TARGETS`). Check that there are no Doxygen
  warnings: CI requires docs on every public member, which the old branch lacked.
- **PR0:** existing test suite unchanged and green; no SoA-path performance regression in md-flexible (LJ, SoA).
- **PR1 non-breaking check:** relative to `origin/master`, PR1 must not change any existing particle, functor or test
  file outside `src/autopas` and `docs`. The whole existing test suite must pass unchanged.
- **New PR1 tests** (`tests/testAutopas/tests/utils/SoATest.cpp` + new files):
  - `sizeof(SoA<Particle::SoAArraysType>) == sizeof(SoAPartition<...>)`.
  - Concept checks: the toy particle satisfies `MultiComponentParticleType`; `ParticleBase`-derived particles don't.
  - Depth growth/shrink; `append` of SoAs and views with different depths; `clear`, `swap`, `pop_back`.
  - `SoAView` component offsets.
  - `SoALoader`/`SoAExtractor` round trip with absent entries: value-initialised entries, and `dummy` ownership.
  - Skip-resize helper with parallel loading.
  - **AoS-vs-SoA equality for the toy functor across all container/traversal/N3 configurations.** Model it on the
    existing Newton3OnOff/AoSvsSoA tests. This covers LCC01, LCC04, VL, VLC and VCL with additional partitions.
- **PR2:**
  - Existing multisite functor tests in `applicationLibrary/molecularDynamics/tests` (AoS vs SoA, N3 on/off,
    mixing on/off).
  - md-flexible multisite example: AoS and SoA energies/trajectories must agree.
  - Timing comparison of the old and new SoA kernels on a multisite example input, recorded in the PR description.
