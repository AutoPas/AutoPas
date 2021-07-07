/**
 * @file OctreeLogger.h
 * @author Johannes Spies
 * @date 21.04.2021
 */

#pragma once

#include <array>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "autopas/containers/octree/OctreeNodeInterface.h"
#include "autopas/containers/octree/OctreeNodeWrapper.h"

namespace autopas {
/**
 * Log an octree to a .vtk file
 * @tparam Particle The enclosed particle type
 */
template <typename Particle>
class OctreeLogger {
 public:
  /**
   * Constructor
   */
  explicit OctreeLogger() = default;

  /**
   * Destructor
   */
  ~OctreeLogger() = default;

  /**
   * Write the octree below the wrapper to a .vtk file
   * @param wrapper A pointer to the octree node wrapper
   */
  void logTree(OctreeNodeWrapper<Particle> *wrapper) { logTree(wrapper->getRaw()); }

  /**
   * This function writes the octree to a .vtk file
   * @param root A pointer to the octree root node
   */
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
    for (int i = 0; i < boxes.size(); ++i) {
      vtkFile << "12\n";  // Write VTK_HEXAHEDRON type for each cell
    }

    // Cleanup
    vtkFile.close();

    // TODO(johannes): Enclose with macro
    //#ifdef AUTOPAS_LOG_OCTREE
    //#endif
  }

  /**
   * Convert a list of octree leaves to JSON and write it to an output file. The output list consists of JSON objects
   * containing the following fields:
   * - `"minmax"`\n
   *   A list of numbers specifying the concatenated minimum and maximum coordinate of the leaf box in 3D.
   *   This is how the coordinate list is encoded in JSON: [x1,y1,z1,x2,y2,z2].
   * - `"fn"`\n
   *   A list of min/max coordinates of all greater than or equal face-neighbors of the leaf
   * - `"fnl"`\n
   *   A list of min/max coordinates of all leaves that touch this leaf via a face.
   * - `"en"`\n
   *   A list of min/max coordinates of all greater than or equal edge-neighbors of the leaf
   * - `"enl"`\n
   *   A list of min/max coordinates of all leaves that touch this leaf via a edge.
   * - `"vn"`\n
   *   A list of min/max coordinates of all greater than or equal vertex-neighbors of the leaf
   * - `"vnl"`\n
   *   A list of min/max coordinates of all leaves that touch this leaf via a vertex.
   *
   * @param out A FILE pointer to the file that should contain the JSON data after the operation
   * @param leaves A list of octree leaves that are echoed into the JSON file
   * @return The FILE pointer is just passed through
   */
  static FILE *leavesToJSON(FILE *out, std::vector<OctreeLeafNode<Particle> *> &leaves) {
    if (out) {
      fprintf(out, "[\n");
      for (int leafIndex = 0; leafIndex < leaves.size(); ++leafIndex) {
        auto leaf = leaves[leafIndex];
        fprintf(out, "{\"minmax\": ");
        outBoxCoordinatesJSON(out, leaf);

        // Print face neighbors
        fprintf(out, ", \"fn\": [");
        bool first = true;
        for (Face *face = getFaces(); *face != O; ++face) {
          auto neighbor = leaf->GTEQ_FACE_NEIGHBOR(*face);
          if (neighbor) {
            if (!first) {
              fprintf(out, ", ");
            }
            first = false;
            outBoxCoordinatesJSON(out, neighbor);
          }
        }

        // Print face neighbor leaves
        fprintf(out, "], \"fnl\": [");
        first = true;
        for (Face *face = getFaces(); *face != O; ++face) {
          auto neighbor = leaf->GTEQ_FACE_NEIGHBOR(*face);
          if (neighbor) {
            auto neighborLeaves = neighbor->getNeighborLeaves(*face);
            for (auto neighborLeaf : neighborLeaves) {
              if (!first) {
                fprintf(out, ", ");
              }
              first = false;
              outBoxCoordinatesJSON(out, neighborLeaf);
            }
          }
        }

        // Print edge neighbors
        fprintf(out, "], \"en\": [");
        first = true;
        for (Edge *edge = getEdges(); *edge != OO; ++edge) {
          auto neighbor = leaf->GTEQ_EDGE_NEIGHBOR(*edge);
          if (neighbor) {
            if (!first) {
              fprintf(out, ", ");
            }
            first = false;
            outBoxCoordinatesJSON(out, neighbor);
          }
        }

        // Print edge neighbor leaves
        fprintf(out, "], \"enl\": [");
        first = true;
        for (Edge *edge = getEdges(); *edge != OO; ++edge) {
          auto neighbor = leaf->GTEQ_EDGE_NEIGHBOR(*edge);
          if (neighbor) {
            auto neighborLeaves = neighbor->getNeighborLeaves(*edge);
            for (auto neighborLeaf : neighborLeaves) {
              if (!first) {
                fprintf(out, ", ");
              }
              first = false;
              outBoxCoordinatesJSON(out, neighborLeaf);
            }
          }
        }

        // Print vertex neighbors
        fprintf(out, "], \"vn\": [");
        first = true;
        for (Vertex *vertex = VERTICES(); *vertex != OOO; ++vertex) {
          auto neighbor = leaf->GTEQ_VERTEX_NEIGHBOR(*vertex);
          if (neighbor) {
            if (!first) {
              fprintf(out, ", ");
            }
            first = false;
            outBoxCoordinatesJSON(out, neighbor);
          }
        }

        // Print vertex neighbor leaves
        fprintf(out, "], \"vnl\": [");
        first = true;
        for (Vertex *vertex = VERTICES(); *vertex != OOO; ++vertex) {
          auto neighbor = leaf->GTEQ_VERTEX_NEIGHBOR(*vertex);
          if (neighbor) {
            auto neighborLeaves = neighbor->getNeighborLeaves(*vertex);

            for (auto neighborLeaf : neighborLeaves) {
              if (!first) {
                fprintf(out, ", ");
              }
              first = false;
              outBoxCoordinatesJSON(out, neighborLeaf);
            }
          }
        }

        fprintf(out, "]}");
        if (leafIndex < (leaves.size() - 1)) {
          fprintf(out, ",");
        }
        fprintf(out, "\n");
      }
      fprintf(out, "]\n");
    } else {
      fprintf(stderr, "ERROR: Dump file is nullptr.\n");
    }

    return out;
  }

  /**
   * Print a list of particle positions to a file as JSON. The list is obtained from the octree root node.
   * @param out The FILE pointer to write the JSON to
   * @param jsonFieldPrefix A short prefix, that will be prepended to the JSON array
   * @param root The root from which the particles should be obtained
   * @return The file pointer
   */
  static FILE *particlesToJSON(FILE *out, OctreeNodeInterface<Particle> *root) {
    if (out) {
      // Get all leaves
      std::vector<OctreeLeafNode<Particle> *> leaves;
      root->appendAllLeaves(leaves);

      fprintf(out, "[");
      for (int leafIndex = 0; leafIndex < leaves.size(); ++leafIndex) {
        OctreeLeafNode<Particle> *leaf = leaves[leafIndex];

        auto n = leaf->numParticles();
        for (int particleIndex = 0; particleIndex < n; ++particleIndex) {
          Particle &particle = leaf->at(particleIndex);
          auto p = particle.getR();
          fputc('[', out);
          out3(out, p[0], p[1], p[2]);
          fputc(']', out);
          if ((particleIndex < (n - 1)) || (leafIndex < (leaves.size() - 1))) {
            fprintf(out, ",");
          }
        }
      }
      fprintf(out, "]");
    } else {
      fprintf(stderr, "ERROR: Dump file is nullptr.\n");
    }

    return out;
  }

  static void octreeToJSON(OctreeNodeInterface<Particle> *owned, OctreeNodeInterface<Particle> *halo,
                           std::vector<OctreeLeafNode<Particle> *> &ownedLeaves,
                           std::vector<OctreeLeafNode<Particle> *> &haloLeaves) {
#if AUTOPAS_LOG_OCTREE
    // Log all owned leaves for this octree
    fclose(OctreeLogger::leavesToJSON(fopen("owned.json", "w"), ownedLeaves));
    // Log all halo leaves for this octree
    fclose(OctreeLogger::leavesToJSON(fopen("halo.json", "w"), haloLeaves));
    FILE *particles = fopen("particles.json", "w");
    fprintf(particles, "{\"owned\": ");
    OctreeLogger::particlesToJSON(particles, owned);
    fprintf(particles, ", \"halo\": \n");
    OctreeLogger::particlesToJSON(particles, halo);
    fprintf(particles, "}");
    fclose(particles);
#endif
  }

 private:
  /**
   *
   * @param out
   * @param x
   * @param y
   * @param z
   */
  static void out3(FILE *out, double x, double y, double z) { fprintf(out, "%.3f,%.3f,%.3f", x, y, z); }

  /**
   * Print the box minimum and maximum coordinates to a given FILE pointer as a JSON list of the form
   * `[min_x, min_y, min_z, max_x, max_y, max_z]`.
   * @param out The FILE pointer
   * @param min The minimum coordinate
   * @param max The maximum coordinate
   */
  static void outLocationArrayJSON(FILE *out, std::array<double, 3> min, std::array<double, 3> max) {
    fputc('[', out);
    out3(out, min[0], min[1], min[2]);
    fputc(',', out);
    out3(out, max[0], max[1], max[2]);
    fputc(']', out);
  }

  /**
   * Print the box minimum and maximum coordinates of an octree node.
   * @param out The FILE pointer
   * @param node An octree node to obtain the box minimum and maximum coordinates from
   */
  static void outBoxCoordinatesJSON(FILE *out, OctreeNodeInterface<Particle> *node) {
    std::array<double, 3> min = node->getBoxMin();
    std::array<double, 3> max = node->getBoxMax();
    outLocationArrayJSON(out, min, max);
  }

  /**
   * Count the iterations to give the written octrees unique filenames
   */
  int unsigned iteration = 0;
};
}  // namespace autopas
