///**
// TODO: These tests are imho not needed anymore. Please confirm in review.
// * @file ParticleIteratorTest.h
// * @author F. Gratl
// * @date 08.03.21
// */
//
// #pragma once
//
// #include <gtest/gtest.h>
//
// #include "AutoPasTestBase.h"
//
// class ParticleIteratorTest : public AutoPasTestBase, public ::testing::WithParamInterface<std::tuple<size_t, size_t>>
// {
// public:
//  /**
//   * Converter functor for generated test names.
//   */
//  struct PrintToStringParamName {
//    template <class ParamType>
//    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
//      auto [numThreads, numAdditionalParticleVectors] = static_cast<ParamType>(info.param);
//      std::stringstream ss;
//      ss << "Threads" << numThreads << "_AdditionalParticleVectors" << numAdditionalParticleVectors;
//      return ss.str();
//    }
//  };
//
// protected:
//  /**
//   * Checks that all particles in a vector of cells are accessed by the iterator.
//   * @param cellsWithParticles Indices of cells that shall contain particles.
//   * @param numAdditionalParticleVectors Number of additional particle vectors the (parallel) loop should consider
//   * besides the cells
//   * @param numThreads Number of threads the test should use.
//   */
//  void testAllParticlesFoundPattern(const std::vector<size_t> &cellsWithParticles, size_t numThreads,
//                                    size_t numAdditionalParticleVectors);
//
//  /**
//   * Applies two iterators and checks that they are behaving exactly the same.
//   * @tparam Iter
//   * @param iter1
//   * @param iter2
//   * @return
//   */
//  template <class Iter>
//  auto iteratorsBehaveEqually(Iter &iter1, Iter &iter2);
//};
