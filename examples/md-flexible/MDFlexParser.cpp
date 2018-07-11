/**
 * @file MDFlexParser.cpp
 * @date 23.02.2018
 * @author F. Gratl
 */

#include "MDFlexParser.h"

bool MDFlexParser::parseInput(int argc, char **argv) {
  bool displayHelp = false;
  int option, option_index;
  constexpr size_t valueOffset = 20;
  static struct option long_options[] = {
      {"container", required_argument, nullptr, 'c'},
      {"cutoff", required_argument, nullptr, 'C'},
      {"data-layout", required_argument, nullptr, 'd'},
      {"functor", required_argument, nullptr, 'f'},
      {"iterations", required_argument, nullptr, 'i'},
      {"particles-per-dimension", required_argument, nullptr, 'n'},
      {"particle-spacing", required_argument, nullptr, 's'},
      {"traversal", required_argument, nullptr, 't'},
  };
  int numOptions = sizeof(long_options) / sizeof(long_options[0]) * 2 + 1;
  if (argc == numOptions) {
    string strArg;
    while ((option = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
      strArg = optarg;
      transform(strArg.begin(), strArg.end(), strArg.begin(), ::tolower);
      switch (option) {
        case 'c': {
          cout << setw(valueOffset) << left << "Container"
               << ":  ";
          if (strArg.find("direct") != string::npos) {
            cout << "DirectSum" << endl;
            containerOption = autopas::directSum;
          } else if (strArg.find("linked") != string::npos) {
            cout << "LinkedCells" << endl;
            containerOption = autopas::linkedCells;
          } else if (strArg.find("verlet") != string::npos) {
            cout << "VerletLists" << endl;
            containerOption = autopas::verletLists;
          } else {
            cerr << "Unknown container : " << strArg << endl;
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
            break;
          }
          cout << setw(valueOffset) << left << "Cutoff radius"
               << ":  " << cutoff << endl;
          break;
        }
        case 'd': {
          cout << setw(valueOffset) << left << "Data Layout"
               << ":  ";
          if (strArg.find("aos") != string::npos) {
            cout << "Array-of-Structures" << endl;
            dataLayoutOption = autopas::aos;
          } else if (strArg.find("soa") != string::npos) {
            cout << "Structure-of-Arrays" << endl;
            dataLayoutOption = autopas::soa;
          } else {
            cerr << "Unknown data layout : " << strArg << endl;
            cerr << "Please use 'AoS' or 'SoA'!" << endl;
            displayHelp = true;
          }
          break;
        }
        case 'f': {
          cout << setw(valueOffset) << left << "Functor"
               << ":  ";
          if (strArg.find("lj") != string::npos || strArg.find("lennard-jones") != string::npos) {
            cout << "Lennard-Jones (12-6)" << endl;
            functorOption = lj12_6;
          } else {
            cerr << "Unknown functor : " << strArg << endl;
            cerr << "Please use 'Lennard-Jones', you have no options here :P" << endl;
            displayHelp = true;
          }
          break;
        }
        case 'i': {
          try {
            iterations = stoul(strArg);
          } catch (const exception &) {
            cerr << "Error parsing number of iterations: " << optarg << endl;
            displayHelp = true;
            break;
          }
          cout << setw(valueOffset) << left << "Iterations"
               << ":  " << iterations << endl;
          break;
        }
        case 'n': {
          try {
            particlesPerDim = stoul(strArg);
          } catch (const exception &) {
            cerr << "Error parsing number of particles per dimension: " << optarg << endl;
            displayHelp = true;
            break;
          }
          cout << "Particles" << endl;
          cout << setw(valueOffset) << left << "  per dimension"
               << ":  " << particlesPerDim << endl;
          cout << setw(valueOffset) << left << "  total"
               << ":  " << (particlesPerDim * particlesPerDim * particlesPerDim) << endl;
          break;
        }
        case 's': {
          try {
            particleSpacing = stod(strArg);
          } catch (const exception &) {
            cerr << "Error parsing separation of particles: " << optarg << endl;
            displayHelp = true;
            break;
          }
          cout << setw(valueOffset) << left << "Particle spacing"
               << ":  " << particleSpacing << endl;
          break;
        }
        case 't': {
          cout << setw(valueOffset) << left << "Allowed traversals"
               << ":  ";
          if (strArg.find("c08") != string::npos) {
            cout << "c08, ";
            traversalOptions.push_back(autopas::TraversalOptions::c08);
          }
          if (strArg.find("slice") != string::npos) {
            cout << "sliced, ";
            traversalOptions.push_back(autopas::TraversalOptions::sliced);
          }
          if (traversalOptions.empty()) {
            cerr << "Unknown Traversal : " << strArg << endl;
            cerr << "Please use 'c08' or 'sliced'!" << endl;
            displayHelp = true;
          }
          cout << "\b\b  " << endl;
          break;
        }
        default: {
          // error message handled by getopt
          displayHelp = true;
        }
      }
    }
  } else {
    cerr << "Wrong number of arguments!" << endl;
    cerr << "Received: " << argc << " Expected: " << numOptions << endl;
    displayHelp = true;
  }
  if (displayHelp) {
    cout << "Usage: " << argv[0] << endl;
    for (auto o : long_options) {
      cout << "    --" << setw(25) << left << o.name;
      if (o.has_arg) cout << "option";
      cout << endl;
    }
    return false;
  }
  return true;
}

autopas::ContainerOptions MDFlexParser::getContainerOption() const { return containerOption; }
double MDFlexParser::getCutoff() const { return cutoff; }
autopas::DataLayoutOption MDFlexParser::getDataLayoutOption() const { return dataLayoutOption; }
MDFlexParser::FunctorOption MDFlexParser::getFunctorOption() const { return functorOption; }

size_t MDFlexParser::getIterations() const { return iterations; }

size_t MDFlexParser::getParticlesPerDim() const { return particlesPerDim; }

double MDFlexParser::getParticleSpacing() const { return particleSpacing; }

const vector<autopas::TraversalOptions> &MDFlexParser::getTraversalOptions() const { return traversalOptions; }
