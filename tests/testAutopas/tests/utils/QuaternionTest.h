/**
* @file QuaternionTest.cpp
* @author S. Newcome
* @date 16/08/2022
*/
#pragma once

#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/Quaternion.h"

/**
* Returns normalized quaternion from direction r and angle theta.
* @param r
* @param theta
* @return
*/
std::array<double, 4> returnNormalizedQuaternion(std::array<double, 3> r, double theta);

std::array<double, 3> returnRotationInAxes(std::array<double,3> pos, int unrotatedAxis, double theta);
