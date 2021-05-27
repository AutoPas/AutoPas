/**
 * @file DomainTools.h
 * @author J. Körner
 * @date 13.05.2021
 */
#pragma once

#include <vector>

namespace DomainTools {
/**
 * Checks if the provided coordinates are inside the domain defined by boxMin and boxMax.
 * @param coordinates The coordinates in question
 * @param boxMin The minimum boundaries of the box domain
 * @param boxMax the maximum boundaries of the box domain
 */
bool isInsideDomain(const std::vector<double> &coordinates, std::vector<double> &boxMin, std::vector<double> &boxMax);

/**
 * Generates a decomposition with a secific number of subdomains.
 * This function uses prime factorization to determine the number of subdomains in each dimension.
 * @param subdomainCount The number of subdomains in the resulting decomposition
 * @param dimensionCount The number of dimensions in the simulation.
 * @param oDecomposition The resulting decomposition.
 */
void generateDecomposition(unsigned int subdomainCount, int dimensionCount, std::vector<double> &oDecomposition);
}
