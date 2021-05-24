/**
 * @file SingleDomain.h
 * @author J. Körner
 * @date 06.05.2021
 */
#pragma once

#include "DomainDecomposition.h"

class SingleDomain final : public DomainDecomposition {
 public:
  SingleDomain(int argc, char **argv, const int &dimensionCount, const std::vector<double> &globalBoxMin,
               const std::vector<double> &globalBoxMax);
  ~SingleDomain() = default;

  void update() override;
  const int getDimensionCount() override { return _dimensionCount; }
  std::vector<double> getGlobalBoxMin() override { return _globalBoxMin; }
  std::vector<double> getGlobalBoxMax() override { return _globalBoxMax; }
  std::vector<double> getLocalBoxMin() override { return _globalBoxMin; }
  std::vector<double> getLocalBoxMax() override { return _globalBoxMax; }
  bool isInsideLocalDomain(const std::vector<double> &coordinates) override;

 private:
  int _dimensionCount;
  std::vector<double> _globalBoxMin;
  std::vector<double> _globalBoxMax;

  void initializeGlobalBox(const std::vector<double> &globalBoxMin, const std::vector<double> &globalBoxMax);
};
