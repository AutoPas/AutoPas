/**
 * @file RDF.h
 * @author D. Martin
 * @date 07.03.2025
 */

#pragma once

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace mdLib {
/**
 * Stores potential values and derivatives to be used in the LuTFunctor
 */
class LookupTable {
 public:
  LookupTable(const std::vector<std::pair<double, double>> &data) : potentialValues(data) {
    // Sort the table according to the distances (r-values), if not already sorted
    std::sort(potentialValues.begin(), potentialValues.end());
  }

  LookupTable() {}

  // Linear interpolation for a given distance r
  double interpolate(double r) const {
    if (potentialValues.empty()) {
      throw std::runtime_error("Lookup table is empty!");
    }

    // If r is outside the range, use the boundary values
    if (r <= potentialValues.front().first) return potentialValues.front().second;
    if (r >= potentialValues.back().first) return potentialValues.back().second;

    // Find the neighboring points
    auto upper = std::lower_bound(
        potentialValues.begin(), potentialValues.end(), std::make_pair(r, 0.0),
        [](const std::pair<double, double> &a, const std::pair<double, double> &b) { return a.first < b.first; });
    auto lower = upper - 1;

    // Linear interpolation
    double r1 = lower->first, U1 = lower->second;
    double r2 = upper->first, U2 = upper->second;
    return U1 + (U2 - U1) * (r - r1) / (r2 - r1);
  }

  // Update the lookup table at runtime
  void updateTable(const std::vector<std::pair<double, double>> &newData) {
    potentialValues = newData;
    std::sort(potentialValues.begin(), potentialValues.end());
  }

  double &operator[](size_t index) {
    if (index >= potentialValues.size()) {
      throw std::out_of_range("LookupTable index out of range!");
    }
    return potentialValues[index].second;
  }

  void writeToCSV(std::string outputFolder, std::string filename) {
    // Ensure the output folder exists
    std::filesystem::create_directories(outputFolder);

    // Construct the full file path
    std::string filePath = outputFolder + "/" + filename + ".csv";

    std::ofstream file(filePath);
    if (!file.is_open()) {
      std::cerr << "Error: Could not open file " << filePath << " for writing." << std::endl;
      return;
    }

    // Write header
    file << "distance,value\n";

    // Write data
    for (const auto &[distance, value] : potentialValues) {
      file << distance << "," << value << "\n";
    }

    file.close();
  }

  void loadFromCSV(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Failed to open file for reading: " + filename);
    }

    std::vector<std::pair<double, double>> newData;
    std::string line;
    while (std::getline(file, line)) {
      std::stringstream ss(line);
      double r, U;
      char comma;
      if (ss >> r >> comma >> U && comma == ',') {
        newData.emplace_back(r, U);
      }
    }
    file.close();
    updateTable(newData);
  }

 private:
  std::vector<std::pair<double, double>> potentialValues;
  std::vector<std::pair<double, double>> derivativeValues;
};
}  // namespace mdLib