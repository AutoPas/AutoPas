/**
 * @file DomainTools.h
 * @author J. Körner
 * @date 13.05.2021
 */
#pragma once

#include <array>

/**
 * Provides some functionality commonly used around rectangular domains.
 */
namespace DomainTools {
/**
 * Checks if the provided coordinates are inside the domain defined by boxMin and boxMax.
 * @param coordinates The coordinates in question.
 * @param boxMin The minimum boundaries of the box domain.
 * @param boxMax the maximum boundaries of the box domain.
 * @return true if the coordinates lie within the provided box.
 */
bool isInsideDomain(const std::array<double, 3> &coordinates, std::array<double, 3> &boxMin,
                    std::array<double, 3> &boxMax);

/**
 * Calculates the distance of provided coordinates to the domain defined by boxMin and boxMax.
 * If the coordinates lie inside the domain, the function returns 0.
 * @param coordinates The coordinates in question.
 * @param boxMin The minimum boundaries of the box domain.
 * @param boxMax the maximum boundaries of the box domain.
 * @return the distance of the coordinates to the provided  domain.
 */
double getDistanceToDomain(const std::array<double, 3> &coordinates, std::array<double, 3> &boxMin,
                           std::array<double, 3> &boxMax);

/**
 * Generates a decomposition with a specific number of subdomains.
 * This function uses prime factorization to determine the number of subdomains in each dimension.
 * @param subdomainCount The number of subdomains in the resulting decomposition.
 * @param decomposition Array containing the number of subdomains per dimension.
 */
void generateDecomposition(unsigned int subdomainCount, std::array<int, 3> &decomposition);

/**
 * Balances two domains by shifting their shared boundary.
 * The balancing depends on work performed in each domain and on the areas of the domain.
 * The two domains have to be part of a regular grid and adjacent to each other.
 * @param leftDomainsWork: The work performed in the domain on the left side of the boundary.
 * @param rightDomainsWork: The work performed in the domain on the right side of the boundary.
 * @param leftDomainsMinBoundaryPosition: The the position of the left domain minimum boundary.
 * @param rightDomainsMaxBoundaryPosition: The the position of the right domain maximum boundary.
 * @returns the updated position of the shared boundary between the two domains.
 */
double balanceAdjacentDomains(const double &leftDomainsWork, const double &rightDomainsWork,
                              const double &leftDomainsMinBoundaryPosition,
                              const double &rightDomainsMaxBoundaryPosition);
}  // namespace DomainTools
