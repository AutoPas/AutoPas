/**
 * @file RDF.h
 * @author D. Martin
 * @date 07.03.2025
 */

#pragma once

#include <algorithm>
#include <cmath>
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
  double interpolateLinearInR(double rSquared) const {
    if (potentialValues.empty()) {
      throw std::runtime_error("Lookup table is empty!");
    }

    double r = std::sqrt(rSquared);

    // Find neighbors: upper is the first point with sqrt(r²) >= r
    auto upper =
        std::lower_bound(potentialValues.begin(), potentialValues.end(), r,
                         [](const std::pair<double, double> &a, double rVal) { return std::sqrt(a.first) < rVal; });

    if (upper == potentialValues.begin()) {
      return upper->second;  // left of range
    }

    if (upper == potentialValues.end()) {
      return (potentialValues.end() - 1)->second;  // right of range
    }

    auto lower = upper - 1;

    double r1 = std::sqrt(lower->first), U1 = lower->second;
    double r2 = std::sqrt(upper->first), U2 = upper->second;

    return U1 + (U2 - U1) * (r - r1) / (r2 - r1);
  }

  // Update the lookup table at runtime
  void updateTable(const std::vector<std::pair<double, double>> &newData) {
    potentialValues = newData;
    std::sort(potentialValues.begin(), potentialValues.end());
  }

  // Computes dU/dr from U(r²), storing the result as (r², dU/dr)
  void computeDerivatives() {
    derivativeValues.clear();
    if (potentialValues.size() < 2) return;

    for (size_t i = 0; i < potentialValues.size() - 1; ++i) {
      double r1Squared = potentialValues[i].first;
      double U1 = potentialValues[i].second;
      double r2Squared = potentialValues[i + 1].first;
      double U2 = potentialValues[i + 1].second;

      // we store derivatives as dUdr2
      double slope = (U2 - U1) / (r2Squared - r1Squared);

      // Store derivative at left endpoint (r1Squared)
      derivativeValues.emplace_back(r1Squared, slope);
    }

    // replicate last slope for the last point
    derivativeValues.emplace_back(potentialValues.back().first, derivativeValues.back().second);
  }

  double getDerivativeValue(double rSquared) const {
    if (derivativeValues.empty()) {
      throw std::runtime_error("Derivative table is empty!");
    }

    // Clamp to bounds
    if (rSquared <= derivativeValues.front().first) return derivativeValues.front().second;
    if (rSquared >= derivativeValues.back().first) return derivativeValues.back().second;

    auto upper = std::lower_bound(
        derivativeValues.begin(), derivativeValues.end(), std::make_pair(rSquared, 0.0),
        [](const std::pair<double, double> &a, const std::pair<double, double> &b) { return a.first < b.first; });

    auto lower = upper - 1;
    return lower->second;  // constant slope in this interval
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

  void writeDerivativesToCSV(const std::string &outputFolder, const std::string &filename) const {
    // Ensure the output folder exists
    std::filesystem::create_directories(outputFolder);

    // Construct the full file path
    std::string filePath = outputFolder + "/" + filename + "_derivatives.csv";

    std::ofstream file(filePath);
    if (!file.is_open()) {
      std::cerr << "Error: Could not open file " << filePath << " for writing." << std::endl;
      return;
    }

    // Write CSV header
    file << "r_squared,dU_dr_squared\n";

    // Write data
    for (const auto &[rSquared, dUdrSquared] : derivativeValues) {
      file << rSquared << "," << dUdrSquared << "\n";
    }

    file.close();
  }

 private:
  std::vector<std::pair<double, double>> potentialValues;
  std::vector<std::pair<double, double>> derivativeValues;
};
}  // namespace mdLib