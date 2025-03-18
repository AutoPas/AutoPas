#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

class LookupTable {
 public:
  LookupTable(const std::vector<std::pair<double, double>> &data) : table(data) {
    // Sortiere die Tabelle nach den Abständen (r-Werten), falls nicht bereits sortiert
    std::sort(table.begin(), table.end());
  }

  LookupTable() {}

  // Lineare Interpolation für einen gegebenen Abstand r
  double interpolate(double r) const {
    if (table.empty()) {
      throw std::runtime_error("Lookup table is empty!");
    }

    // Falls r außerhalb des Bereichs liegt, nutze die Randwerte
    if (r <= table.front().first) return table.front().second;
    if (r >= table.back().first) return table.back().second;

    // Finde die benachbarten Punkte
    auto upper = std::lower_bound(
        table.begin(), table.end(), std::make_pair(r, 0.0),
        [](const std::pair<double, double> &a, const std::pair<double, double> &b) { return a.first < b.first; });
    auto lower = upper - 1;

    // Lineare Interpolation
    double r1 = lower->first, U1 = lower->second;
    double r2 = upper->first, U2 = upper->second;
    return U1 + (U2 - U1) * (r - r1) / (r2 - r1);
  }

  // Update der Lookup-Tabelle zur Laufzeit
  void updateTable(const std::vector<std::pair<double, double>> &newData) {
    table = newData;
    std::sort(table.begin(), table.end());
  }

  double &operator[](size_t index) {
    if (index >= table.size()) {
      throw std::out_of_range("LookupTable index out of range!");
    }
    return table[index].second;
  }

  void saveToCSV(const std::string &filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Failed to open file for writing: " + filename);
    }
    for (const auto &entry : table) {
      file << entry.first << "," << entry.second << "\n";
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
  std::vector<std::pair<double, double>> table;
};