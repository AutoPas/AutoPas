/**
 * @file StaticSelectorMacros.h
 * @author seckler
 * @date 21.06.18
 */

#pragma once

// either LinkedCells<Particle,ParticleCell>
// or VerletLists<Particle>
// or DirectSum<Particle, ParticleCell>
/**
 * Will execute the passed body with the static container type of container, i.e. either
 * LinkedCells, VerletLists or DirectSum
 * @param container the container to be used
 * @param body the function body to be executed
 */
#define WithStaticContainerType(container, body)                                                                    \
  {                                                                                                                 \
    auto container_ptr = container.get();                                                                           \
    if (auto container = dynamic_cast<                                                                              \
            autopas::LinkedCells<typename std::remove_pointer_t<decltype(container_ptr)>::ParticleType,             \
                                 typename std::remove_pointer_t<decltype(container_ptr)>::ParticleCellType>*>(      \
            container_ptr)) {                                                                                       \
      body                                                                                                          \
    } else if (auto container = dynamic_cast<                                                                       \
                   autopas::VerletLists<typename std::remove_pointer_t<decltype(container_ptr)>::ParticleType>*>(   \
                   container_ptr)) {                                                                                \
      body                                                                                                          \
    } else if (auto container = dynamic_cast<                                                                       \
                   autopas::DirectSum<typename std::remove_pointer_t<decltype(container_ptr)>::ParticleType,        \
                                      typename std::remove_pointer_t<decltype(container_ptr)>::ParticleCellType>*>( \
                   container_ptr)) {                                                                                \
      body                                                                                                          \
    } else {                                                                                                        \
      autopas::utils::ExceptionHandler::exception("wrong type of container in StaticSelectorMacros.h");             \
    }                                                                                                               \
  }
