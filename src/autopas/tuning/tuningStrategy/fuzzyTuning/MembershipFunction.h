/**
 * @file MembershipFunction.h
 * @author Manuel Lerchner
 * @date 17.04.24
 */

#pragma once

#include "FuzzySet.h"

namespace autopas::fuzzy_logic {

/**
 * Constructs a Triangle FuzzySet with the given linguistic term and membership function.
 * @param linguisticTerm The linguistic term of this FuzzySet.
 * @param min The left border of the triangle. The membership function will be 0 for values smaller than min.
 * @param mid The peak of the triangle. The membership function will be 1 for values equal to mid.
 * @param max The right border of the triangle. The membership function will be 0 for values greater than max.
 */
std::shared_ptr<FuzzySet> makeTriangle(const std::string &linguisticTerm, double min, double mid, double max);

/**
 * Constructs a Trapezoid FuzzySet with the given linguistic term and membership function.
 * @param linguisticTerm The linguistic term of this FuzzySet.
 * @param min The left border of the trapezoid. The membership function will be 0 for values smaller than min.
 * @param mid1 The left peak of the trapezoid. The membership function will be 1 for values equal to mid1.
 * @param mid2 The right peak of the trapezoid. The membership function will be 1 for values equal to mid2.
 * @param max The right border of the trapezoid. The membership function will be 0 for values greater than max.
 */
std::shared_ptr<FuzzySet> makeTrapezoid(const std::string &linguisticTerm, double min, double mid1, double mid2,
                                        double max);

/**
 * Constructs a Gaussian FuzzySet with the given linguistic term and membership function.
 * @param linguisticTerm The linguistic term of this FuzzySet.
 * @param mean The mean of the Gaussian. The membership function will be 1 for values equal to mean.
 * @param sigma The standard deviation of the Gaussian.
 */
std::shared_ptr<FuzzySet> makeGaussian(const std::string &linguisticTerm, double mean, double sigma);

}  // namespace autopas::fuzzy_logic