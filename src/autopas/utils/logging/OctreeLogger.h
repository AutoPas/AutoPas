/**
 * @file OctreeLogger.h
 * @author Johannes Spies
 * @date 21.04.2021
 */

#pragma once

#include <string.h>

#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "autopas/containers/octree/OctreeNodeInterface.h"

namespace autopas {
/**
 * Log an octree to a .vtk file
 */
class OctreeLogger {
 public:
  /**
   * Constructor
   */
  explicit OctreeLogger();

  /**
   * Destructor
   */
  ~OctreeLogger();

  /**
   * This function writes the octree to a .vtk file
   * @tparam Particle The enclosed particle type
   * @param root A pointer to the octree root node
   */
  template <typename Particle>
  void logTree(OctreeNodeInterface<Particle> *root) {
    // Load the leaf boxes
    using Position = std::array<double, 3>;
    using Box = std::pair<Position, Position>;
    std::vector<Box> boxes;
    root->appendAllLeafBoxes(boxes);
    auto boxCount = boxes.size();
    auto pointCount = 8 * boxCount;

    // Open the VTK file
    char filename[256] = {0};
    snprintf(filename, sizeof(filename), "octree_%d.vtk", iteration++);
    std::ofstream vtkFile;
    vtkFile.open(filename);

    if (not vtkFile.is_open()) {
      throw std::runtime_error("OctreeLogger::logTree(): Failed to open file \"" + std::string(filename) + "\".");
    }

    // Write the header
    vtkFile << "# vtk DataFile Version 2.0\n"
            << "Octree boxes\n"
            << "ASCII\n"

            << "DATASET UNSTRUCTURED_GRID\n"
            << "\n";

    // Write points
    vtkFile << "POINTS " << pointCount << " float\n";  // Points header
    for (Box box : boxes) {
      Position min = box.first;
      Position max = box.second;

      auto [minX, minY, minZ] = min;
      auto [maxX, maxY, maxZ] = max;

      // Write the points in the order of the VTK_HEXAHEDRON. Each point is
      // written on its own line.

      vtkFile << minX << " " << minY << " " << minZ << "\n";  // 0 ---
      vtkFile << maxX << " " << minY << " " << minZ << "\n";  // 1 +--
      vtkFile << maxX << " " << maxY << " " << minZ << "\n";  // 2 ++-
      vtkFile << minX << " " << maxY << " " << minZ << "\n";  // 3 -+-
      vtkFile << minX << " " << minY << " " << maxZ << "\n";  // 4 --+
      vtkFile << maxX << " " << minY << " " << maxZ << "\n";  // 5 +-+
      vtkFile << maxX << " " << maxY << " " << maxZ << "\n";  // 6 +++
      vtkFile << minX << " " << maxY << " " << maxZ << "\n";  // 7 -++
    }
    vtkFile << "\n";

    // Write cells
    auto cellListSize = pointCount + boxCount;
    vtkFile << "CELLS " << boxCount << " " << cellListSize << "\n";
    for (auto boxIndex = 0; boxIndex < boxCount; ++boxIndex) {
      vtkFile << "8 ";  // Output # of elements in the following line.
      for (auto pointIndex = 0; pointIndex < 8; ++pointIndex) {
        // Generate an index that points to the corresponding point in the points list.
        auto offset = 8 * boxIndex + pointIndex;
        vtkFile << offset;
        if (pointIndex < 7) {
          vtkFile << " ";
        }
      }
      vtkFile << "\n";
    }
    vtkFile << "\n";

    // Write cell types
    vtkFile << "CELL_TYPES " << boxCount << "\n";
    for (Box box : boxes) {
      vtkFile << "12\n";  // Write VTK_HEXAHEDRON type for each cell
    }

    // Cleanup
    vtkFile.close();

    // TODO(johannes): Enclose with macro
    //#ifdef AUTOPAS_LOG_OCTREE
    //#endif
  }

 private:
  /**
   * Count the iterations to give the written octrees unique filenames
   */
  int unsigned iteration = 0;
};
}  // namespace autopas
