#pragma once

namespace ParticleTypes {

constexpr unsigned long FLUID = 0;

constexpr unsigned long WALL_SMOOTH = 1;
constexpr unsigned long WALL_BOSS = 2;
constexpr unsigned long WALL_PIT = 3;
constexpr unsigned long WALL_GRID = 4;
constexpr unsigned long WALL_DUAL_BOSS = 5;

// Wettability parameter zeta (paper model: effective_epsilon_fw = zeta * epsilon_ff).
// Also set ZETA in surface_generator.py to the same value for consistent documentation.
constexpr double ZETA = 0.6;

template <class TypeId>
constexpr bool isWall(TypeId typeId) {
  return typeId == static_cast<TypeId>(WALL_SMOOTH) or typeId == static_cast<TypeId>(WALL_BOSS) or
         typeId == static_cast<TypeId>(WALL_PIT) or typeId == static_cast<TypeId>(WALL_GRID) or
         typeId == static_cast<TypeId>(WALL_DUAL_BOSS);
}

template <class TypeId>
constexpr bool isFluid(TypeId typeId) {
  return typeId == static_cast<TypeId>(FLUID);
}

template <class TypeId>
constexpr bool isFluidWallPair(TypeId firstTypeId, TypeId secondTypeId) {
  return (isFluid(firstTypeId) and isWall(secondTypeId)) or (isWall(firstTypeId) and isFluid(secondTypeId));
}

template <class TypeId>
constexpr bool isWallWallPair(TypeId firstTypeId, TypeId secondTypeId) {
  return isWall(firstTypeId) and isWall(secondTypeId);
}

}  // namespace ParticleTypes
