/**
 * @file MDFlexParser.cpp
 * @date 23.02.2018
 * @author F. Gratl
 */

#include "MDFlexParser.h"

bool MDFlexParser::parseInput(int argc, char **argv) {
  bool displayHelp = false;
  int option, option_index;
  static struct option long_options[] = {{"box-length", required_argument, nullptr, 'b'},
                                         {"container", required_argument, nullptr, 'c'},
                                         {"cutoff", required_argument, nullptr, 'C'},
                                         {"distribution-mean", required_argument, nullptr, 'm'},
                                         {"distribution-stddeviation", required_argument, nullptr, 'z'},
                                         {"data-layout", required_argument, nullptr, 'd'},
                                         {"functor", required_argument, nullptr, 'f'},
                                         {"help", no_argument, nullptr, 'h'},
                                         {"iterations", required_argument, nullptr, 'i'},
                                         {"no-flops", no_argument, nullptr, 'F'},
                                         {"particles-generator", required_argument, nullptr, 'g'},
                                         {"particles-per-dimension", required_argument, nullptr, 'n'},
                                         {"particle-spacing", required_argument, nullptr, 's'},
                                         {"traversal", required_argument, nullptr, 't'},
                                         {"tuning-interval", required_argument, nullptr, 'I'},
                                         {"log-level", required_argument, nullptr, 'l'},
                                         {"verlet-rebuild-frequency", required_argument, nullptr, 'v'},
                                         {"verlet-skin-radius", required_argument, nullptr, 'r'},
                                         {"vtk", required_argument, nullptr, 'w'},
                                         {nullptr, 0, nullptr, 0}};  // needed to signal the end of the array
  string strArg;
  int bla = 0;
  while ((option = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
    ++bla;
    if (optarg != nullptr) strArg = optarg;
    transform(strArg.begin(), strArg.end(), strArg.begin(), ::tolower);
    switch (option) {
      case 'b': {
        try {
          boxLength = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing box length: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'c': {
        if (strArg.find("direct") != string::npos) {
          containerOption = autopas::directSum;
        } else if (strArg.find("linked") != string::npos or strArg.find("lc") != string::npos) {
          containerOption = autopas::linkedCells;
        } else if (strArg.find("verlet") != string::npos or strArg.find("vl") != string::npos) {
          containerOption = autopas::verletLists;
        } else {
          cerr << "Unknown container option: " << strArg << endl;
          cerr << "Please use 'DirectSum' or 'LinkedCells'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'C': {
        try {
          cutoff = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing cutoff Radius: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'd': {
        if (strArg.find("aos") != string::npos) {
          dataLayoutOption = autopas::aos;
        } else if (strArg.find("soa") != string::npos) {
          dataLayoutOption = autopas::soa;
        } else {
          cerr << "Unknown data layout : " << strArg << endl;
          cerr << "Please use 'AoS' or 'SoA'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'f': {
        if (strArg.find("lj") != string::npos || strArg.find("lennard-jones") != string::npos) {
          functorOption = lj12_6;
        } else {
          cerr << "Unknown functor : " << strArg << endl;
          cerr << "Please use 'Lennard-Jones', you have no options here :P" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'F': {
        measureFlops = false;
        break;
      }
      case 'g': {
        if (strArg.find("grid") != string::npos) {
          generatorOption = GeneratorOption::grid;
        } else if (strArg.find("gaus") != string::npos) {
          generatorOption = GeneratorOption::gaussian;
        } else {
          cerr << "Unknown generator : " << strArg << endl;
          cerr << "Please use 'Grid' or 'Gaussian'" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'h': {
        displayHelp = true;
        break;
      }
      case 'i': {
        try {
          iterations = stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing number of iterations: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'I': {
        try {
          tuningInterval = (unsigned int)stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing tuning interval: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'l': {
        switch (strArg[0]) {
          case 't': {
            logLevel = spdlog::level::trace;
            break;
          }
          case 'd': {
            logLevel = spdlog::level::debug;
            break;
          }
          case 'i': {
            logLevel = spdlog::level::info;
            break;
          }
          case 'w': {
            logLevel = spdlog::level::warn;
            break;
          }
          case 'e': {
            logLevel = spdlog::level::err;
            break;
          }
          case 'c': {
            logLevel = spdlog::level::critical;
            break;
          }
          case 'o': {
            logLevel = spdlog::level::off;
            break;
          }
          default: {
            cerr << "Unknown Log Level : " << strArg << endl;
            cerr << "Please use 'trace', 'debug', 'info', 'warning', 'error', 'critical' or 'off'." << endl;
            displayHelp = true;
          }
        }
        break;
      }
      case 'm': {
        try {
          distributionMean = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing distribution mean: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'n': {
        try {
          particlesPerDim = stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing number of particles per dimension: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 's': {
        try {
          particleSpacing = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing separation of particles: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 't': {
        if (strArg.find("c08") != string::npos) {
          traversalOptions.push_back(autopas::TraversalOptions::c08);
        }
        if (strArg.find("sli") != string::npos) {
          traversalOptions.push_back(autopas::TraversalOptions::sliced);
        }
        if (traversalOptions.empty()) {
          cerr << "Unknown Traversal : " << strArg << endl;
          cerr << "Please use 'c08' or 'sliced'!" << endl;
          displayHelp = true;
        }
        break;
      }
      case 'v': {
        try {
          verletRebuildFrequency = stoul(strArg);
        } catch (const exception &) {
          cerr << "Error parsing verlet-rebuild-frequency: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'w': {
        writeVTK = strArg;
        break;
      }
      case 'r': {
        try {
          verletSkinRadius = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing verlet-skin-radius: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      case 'z': {
        try {
          distributionStdDev = stod(strArg);
        } catch (const exception &) {
          cerr << "Error parsing distribution standard deviation: " << optarg << endl;
          displayHelp = true;
        }
        break;
      }
      default: {
        // error message handled by getopt
        displayHelp = true;
      }
    }
  }
  if (displayHelp) {
    cout << "Usage: " << argv[0] << endl;
    for (auto o : long_options) {
      if (o.name == nullptr) continue;
      cout << "    --" << setw(valueOffset + 2) << left << o.name;
      if (o.has_arg) cout << "option";
      cout << endl;
    }
    return false;
  }
  return true;
}

void MDFlexParser::printConfig() {
  constexpr size_t valueOffset = 32;
  cout << setw(valueOffset) << left << "Container"
       << ":  ";
  switch (containerOption) {
    case autopas::ContainerOptions::directSum: {
      cout << "DirectSum" << endl;
      break;
    }
    case autopas::ContainerOptions::linkedCells: {
      cout << "LinkedCells" << endl;
      break;
    }
    case autopas::ContainerOptions::verletLists: {
      cout << "VerletLists" << endl;
      break;
    }
  }

  if (containerOption == autopas::ContainerOptions::verletLists) {
    cout << setw(valueOffset) << left << "Verlet rebuild frequency"
         << ":  " << verletRebuildFrequency << endl;

    cout << setw(valueOffset) << left << "Verlet skin radius"
         << ":  " << verletSkinRadius << endl;
  }

  cout << setw(valueOffset) << left << "Data Layout"
       << ":  ";
  switch (dataLayoutOption) {
    case autopas::DataLayoutOption::aos: {
      cout << "Array-of-Structures" << endl;
      break;
    }
    case autopas::DataLayoutOption::soa: {
      cout << "Structure-of-Arrays" << endl;
      break;
    }
  }

  cout << setw(valueOffset) << left << "Functor"
       << ":  ";
  switch (functorOption) {
    case FunctorOption::lj12_6: {
      cout << "Lennard-Jones (12-6)" << endl;
      break;
    }
  }

  cout << setw(valueOffset) << left << "Cutoff radius"
       << ":  " << cutoff << endl;

  cout << setw(valueOffset) << left << "Particle Generator"
       << ":  ";
  switch (generatorOption) {
    case GeneratorOption::grid: {
      cout << "Grid generator" << endl;
      cout << setw(valueOffset) << left << "Particle spacing"
           << ":  " << particleSpacing << endl;
      break;
    }
    case GeneratorOption::gaussian: {
      cout << "Gaussian generator" << endl;
      cout << setw(valueOffset) << left << "Box length"
           << ":  " << boxLength << endl;
      cout << setw(valueOffset) << left << "Distribution mean"
           << ":  " << distributionMean << endl;
      cout << setw(valueOffset) << left << "Distribution standard deviation"
           << ":  " << distributionStdDev << endl;
      break;
    }
  }

  cout << "Particles" << endl;
  cout << setw(valueOffset) << left << "  per dimension"
       << ":  " << particlesPerDim;
  if (generatorOption != grid) cout << " (approximately)";
  cout << endl;
  cout << setw(valueOffset) << left << "  total"
       << ":  " << (particlesPerDim * particlesPerDim * particlesPerDim) << endl;

  cout << setw(valueOffset) << left << "Allowed traversals"
       << ":  ";
  for (auto &t : traversalOptions) {
    switch (t) {
      case autopas::TraversalOptions::c08: {
        cout << "c08, ";
        break;
      }
      case autopas::TraversalOptions::sliced: {
        cout << "sliced, ";
        break;
      }
    }
  }
  // deletes last comma
  cout << "\b\b  " << endl;

  cout << setw(valueOffset) << left << "Iterations"
       << ":  " << iterations << endl;
  cout << setw(valueOffset) << left << "Tuning Interval"
       << ":  " << tuningInterval << endl;
}

autopas::ContainerOptions MDFlexParser::getContainerOption() const { return containerOption; }

double MDFlexParser::getCutoff() const { return cutoff; }

autopas::DataLayoutOption MDFlexParser::getDataLayoutOption() const { return dataLayoutOption; }

MDFlexParser::FunctorOption MDFlexParser::getFunctorOption() const { return functorOption; }

size_t MDFlexParser::getIterations() const { return iterations; }

size_t MDFlexParser::getParticlesPerDim() const { return particlesPerDim; }

double MDFlexParser::getParticleSpacing() const { return particleSpacing; }

const vector<autopas::TraversalOptions> &MDFlexParser::getTraversalOptions() const { return traversalOptions; }

size_t MDFlexParser::getVerletRebuildFrequency() const { return verletRebuildFrequency; }

double MDFlexParser::getVerletSkinRadius() const { return verletSkinRadius; }

MDFlexParser::GeneratorOption MDFlexParser::getGeneratorOption() const { return generatorOption; }

double MDFlexParser::getDistributionMean() const { return distributionMean; }

double MDFlexParser::getDistributionStdDev() const { return distributionStdDev; }

string MDFlexParser::getWriteVTK() const { return writeVTK; }

double MDFlexParser::getBoxLength() const {
  if (boxLength == -1) return ceil(2 * distributionMean);
  return boxLength;
}

bool MDFlexParser::getMeasureFlops() const { return measureFlops; }

unsigned int MDFlexParser::getTuningInterval() const { return tuningInterval; }

autopas::Logger::LogLevel MDFlexParser::getLogLevel() const { return logLevel; }
//spdlog::level::level_enum MDFlexParser::getLogLevel() const { return logLevel; }
