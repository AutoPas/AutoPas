/**
 * @file DomainTools.h
 * @author J. Körner
 * @date 13.05.2021
 */
#pragma once

#include <array>
#include <cstddef>

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
bool isInsideDomain(const std::array<double, 3> &coordinates, const std::array<double, 3> &boxMin,
                    const std::array<double, 3> &boxMax);

/**
 * Generates a decomposition with a specific number of subdomains.
 * This function uses prime factorization to determine the number of subdomains in each dimension.
 * @param subdomainCount The number of subdomains in the resulting decomposition.
 * @param decomposition Array containing the number of subdomains per dimension.
 */
void generateDecomposition(unsigned int subdomainCount, std::array<int, 3> &decomposition);

/**
 * Converts a domain id to the domain index, i.e. rank of the local processor.
 * @param domainId: the domain id to be converted to an index.
 * @param decomposition: The global domain's decomposition.
 * @return the domain index which corresponds to the provided domain id.
 */
int convertIdToIndex(const std::array<int, 3> &domainId, const std::array<int, 3> decomposition);

/**
 * Convert a domain index to the corresponding domain id.
 * @param domainIndex: The index to be converted to a domain id.
 * @param decomposition: The global domain's decomposition.
 * @return the domain id which corresponds to the provided domain index.
 */
std::array<int, 3> convertIndexToId(const int domainIndex, const std::array<int, 3> decomposition);

/**
 * Returns the accumulated tail of a decomposition.
 * The accumulated tail refers to the number of elements in a decomposition's tail.
 * If the decomposition is [4,3,2], the accumulated tail of the subdivision starting after index 0 is 3 * 2 = 6.
 * The accumulated tail of this decomposition starting after index 1 is 2 .
 * the accumulated tail of this decomposition starting after index 2 is 1.
 * @param index: The starting index of the decomposition's tail.
 * @param decomposition: The global domain's decomposition.
 * @return the accumulated tail starting at index.
 */
int getAccumulatedTail(const size_t index, const std::array<int, 3> decomposition);

/**
 * Returns the extent of a subdomain in a domain decomposition.
 * The extent is defined by the starting coordinates of the subdomain within the domain decomposition.
 * A subdomain with id [2, 3, 4] in a subdivision [5, 5, 5] has the corresponding extent [2, 3, 3, 4, 4, 5],  where
 * the values are defined as follows: [ xmin, xmax, ymin, ymax, zmin, zmax ]. The global domain extent of this domain
 * decomposition would be [ 0, 5, 0, 5, 0, 5 ].
 * @param subdomainIndex: The index of the subdomain for which to calculate the domain extent.
 * @param decomposition: The decomposition subdomain's respective global domain.
 * @return the extent of the subdomain with index subdomainIndex.
 */
std::array<int, 6> getExtentOfSubdomain(const int subdomainIndex, const std::array<int, 3> decomposition);
}  // namespace DomainTools
