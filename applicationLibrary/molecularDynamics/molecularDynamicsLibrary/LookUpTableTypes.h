/**
 * @file LookUpTableTypes.h
 * @author J. Hampe
 * @date 29.1.2024
 */

#ifndef AUTOPAS_LOOKUPTABLETYPES_H
#define AUTOPAS_LOOKUPTABLETYPES_H

namespace ForceLookUpTable {

enum IntervalType { evenSpacing };

enum InterpolationType { nextNeighbor, linear };

enum PositionType { absolute, relative };
}  // namespace ForceLookUpTable

#endif  // AUTOPAS_LOOKUPTABLETYPES_H
