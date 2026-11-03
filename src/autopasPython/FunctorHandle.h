/**
 * @file FunctorHandle.h
 * @date 03.11.2026
 */
#pragma once

#include <limits>
#include <string>

#include "autopas/AutoPasDecl.h"
#include "autopas/options/InteractionTypeOption.h"
#include "src/TypeDefinitions.h"

namespace autopasPython {

/**
 * The AutoPas type all of the Python interface operates on.
 */
using AutoPasType = autopas::AutoPas<ParticleType>;

/**
 * Type erased view on a functor.
 *
 * AutoPas functors are CRTP classes and AutoPas::computeInteractions() is a template over them. Python, however, has
 * to choose the functor at runtime. Each implementation of this interface therefore belongs to one explicit
 * instantiation of computeInteractions().
 */
class FunctorHandle {
 public:
  virtual ~FunctorHandle() = default;

  /**
   * Whether this functor is pairwise or triwise.
   * @return interaction type.
   */
  [[nodiscard]] virtual autopas::InteractionTypeOption::Value interactionType() const = 0;

  /**
   * Whether this functor is implemented in Python and therefore restricted to the AoS data layout.
   * @return true if and only if the functor calls into Python.
   */
  [[nodiscard]] virtual bool isPython() const = 0;

  /**
   * Whether the functor offers a Newton 3 implementation.
   * @return true if Newton 3 can be used.
   */
  [[nodiscard]] virtual bool allowsNewton3() const = 0;

  /**
   * Whether the functor offers a non Newton 3 implementation.
   * @return true if Newton 3 can be disabled.
   */
  [[nodiscard]] virtual bool allowsNonNewton3() const = 0;

  /**
   * Name of the functor.
   * @return name
   */
  [[nodiscard]] virtual std::string name() const = 0;

  /**
   * Creates the actual AutoPas functor. A handle can only be bound once, hence only be used by one simulation.
   * @param cutoff The cutoff of the simulation.
   * @param particlePropertiesLibrary The library of the simulation. It has to outlive this object.
   * @param relevantForTuning If the AutoTuner should measure this functor. There is one per interaction type.
   */
  virtual void bind(double cutoff, ParticlePropertiesLibraryType &particlePropertiesLibrary,
                    bool relevantForTuning) = 0;

  /**
   * Applies the functor to all particles in the given container.
   * @note The GIL is not held when this is called.
   * @param autoPas The container.
   * @return true if this was a tuning iteration.
   */
  virtual bool compute(AutoPasType &autoPas) = 0;

  /**
   * Potential energy of the latest compute() call.
   * @return potential energy or NaN if this functor does not provide it.
   */
  [[nodiscard]] virtual double potentialEnergy() const { return std::numeric_limits<double>::quiet_NaN(); }

  /**
   * Virial of the latest compute() call.
   * @return virial or NaN if this functor does not provide it.
   */
  [[nodiscard]] virtual double virial() const { return std::numeric_limits<double>::quiet_NaN(); }
};

}  // namespace autopasPython
