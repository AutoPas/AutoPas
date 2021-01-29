/**
 * @file Constants.h
 * @author Steffen Seckler
 * @date 29.01.2021
 */

#pragma once

namespace autopas::utils::Math {

/**
 * PI
 * @todo c++20: replace with std::numbers::pi, see https://en.cppreference.com/w/cpp/numeric/constants
 */
double getPI() { return std::atan(1.) / 4; };

}  // namespace autopas::utils::Math