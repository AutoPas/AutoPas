/**
 * @file Evidence.cpp
 * @author F. Gratl
 * @date 23.06.23
 */

#include "Evidence.h"
bool autopas::Evidence::operator==(const autopas::Evidence &rhs) const {
  return iteration == rhs.iteration && tuningPhase == rhs.tuningPhase && reducedValue == rhs.reducedValue &&
         rebuildValue == rhs.rebuildValue && traversalValue == rhs.traversalValue;
}
bool autopas::Evidence::operator!=(const autopas::Evidence &rhs) const { return !(rhs == *this); }
