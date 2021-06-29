/**
 * @file DomainTools.h
 * @author J. Körner
 * @date 13.05.2021
 */
#pragma once

#include <array>
#include <vector>

namespace DomainTools {
/**
 * Checks if the provided coordinates are inside the domain defined by boxMin and boxMax.
 * @param coordinates The coordinates in question.
 * @param boxMin The minimum boundaries of the box domain.
 * @param boxMax the maximum boundaries of the box domain.
 * $return true if the coordinates lie within the provided box.
 */
bool isInsideDomain(const std::vector<double> &coordinates, std::vector<double> &boxMin, std::vector<double> &boxMax);

/**
 * Checks if the provided coordinates are inside the domain defined by boxMin and boxMax.
 * Instead of a vector, the parameters are of type std::array<double, 3> to be compatible with AutoPas.
 * @param coordinates The coordinates in question.
 * @param boxMin The minimum boundaries of the box domain.
 * @param boxMax the maximum boundaries of the box domain.
 * $return true if the coordinates lie within the provided box.
 */
bool isInsideDomain(const std::array<double, 3> &coordinates, std::array<double, 3> &boxMin,
                    std::array<double, 3> &boxMax);

/**
 * Calculates the distance of provided coordinates to the domain defined by boxMin and boxMax.
 * If the coordinates lie inside the domain, the function returns 0.
 * If the coordinates and the domain do not have the same number of dimensions, the function returns -1.
 * @param coordinates The coordinates in question.
 * @param boxMin The minimum boundaries of the box domain.
 * @param boxMax the maximum boundaries of the box domain.
 * @return the distance of the coordinates to the provided  domain.
 */
double getDistanceToDomain(const std::vector<double> &coordinates, std::vector<double> &boxMin,
                           std::vector<double> &boxMax);

/**
 * Calculates the distance of provided coordinates to the domain defined by boxMin and boxMax.
 * Instead of a vector, the parameters are of type std::array<double, 3> to be compatible with AutoPas.
 * If the coordinates lie inside the domain, the function returns 0.
 * If the coordinates and the domain do not have the same number of dimensions, the function returns -1.
 * @param coordinates The coordinates in question.
 * @param boxMin The minimum boundaries of the box domain.
 * @param boxMax the maximum boundaries of the box domain.
 * @return the distance of the coordinates to the provided  domain.
 */
double getDistanceToDomain(const std::array<double, 3> &coordinates, std::array<double, 3> &boxMin,
                           std::array<double, 3> &boxMax);

/**
 * Generates a decomposition with a secific number of subdomains.
 * This function uses prime factorization to determine the number of subdomains in each dimension.
 * @param subdomainCount The number of subdomains in the resulting decomposition.
 * @param dimensionCount The number of dimensions in the simulation.
 * @param decomposition Vector containing the number of subdomains per dimension.
 */
void generateDecomposition(unsigned int subdomainCount, int dimensionCount, std::vector<int> &decomposition);
}  // namespace DomainTools
