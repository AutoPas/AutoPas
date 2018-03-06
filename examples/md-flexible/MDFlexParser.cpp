#include "MDFlexParser.h"

bool MDFlexParser::parseInput(int argc, char **argv) {
  bool displayHelp = false;
  int option, option_index;
  static struct option long_options[] = {
      {"container", required_argument, nullptr, 'c'},
      {"cutoff", required_argument, nullptr, 'C'},
      {"data-layout", required_argument, nullptr, 'd'},
      {"functor", required_argument, nullptr, 'f'},
      {"iterations", required_argument, nullptr, 'i'},
      {"particles-per-dimension", required_argument, nullptr, 'n'},
      {"particle-spacing", required_argument, nullptr, 's'},
  };
  int numOptions = sizeof(long_options) / sizeof(long_options[0]) * 2 + 1;
  if (argc == numOptions) {
    string strArg;
    while ((option = getopt_long(argc, argv, "", long_options,
                                 &option_index)) != -1) {
      strArg = optarg;
      transform(strArg.begin(), strArg.end(), strArg.begin(), ::tolower);
      switch (option) {
        case 'c': {
          if (strArg.find("direct") != string::npos) {
            cout << "Using container: DirectSum" << endl;
            containerOption = directSum;
          } else if (strArg.find("linked") != string::npos) {
            cout << "Using container: LinkedCells" << endl;
            containerOption = linkedCells;
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
          cout << "Cutoff radius: " << cutoff << endl;
          break;
        }
        case 'd': {
          if (strArg.find("aos") != string::npos) {
            cout << "Using Array-of-Structures" << endl;
            dataLayoutOption = aos;
          } else if (strArg.find("soa") != string::npos) {
            cout << "Using Structure-of-Arrays" << endl;
            dataLayoutOption = soa;
          } else {
            cerr << "Unknown data layout : " << strArg << endl;
            cerr << "Please use 'AoS' or 'SoA'!" << endl;
            displayHelp = true;
          }
          break;
        }
        case 'f': {
          if (strArg.find("lj") != string::npos ||
              strArg.find("lennard-jones") != string::npos) {
            cout << "Using Lennard-Jones (12-6) Functor" << endl;
            functorOption = lj12_6;
          } else {
            cerr << "Unknown functor : " << strArg << endl;
            cerr << "Please use 'Lennard-Jones', you have no options here :P"
                 << endl;
            displayHelp = true;
          }
          break;
        }
        case 'i': {
          try {
            iterations = stoul(strArg);
          } catch (const exception &) {
            cerr << "Error parsing number of iterations: "
                 << optarg << endl;
            displayHelp = true;
            break;
          }
          cout << "Simulating " << iterations
               << " iterations." << endl;
          break;
        }
        case 'n': {
          try {
            particlesPerDim = stoul(strArg);
          } catch (const exception &) {
            cerr << "Error parsing number of particles per dimension: "
                 << optarg << endl;
            displayHelp = true;
            break;
          }
          cout << "Simulating " << particlesPerDim
               << " particles per dimension." << endl;
          break;
        }
        case 's': {
          try {
            particleSpacing = stod(strArg);
          } catch (const exception &) {
            cerr << "Error parsing separation of particles: "
                 << optarg << endl;
            displayHelp = true;
            break;
          }
          cout << "Particles are separated by " << particleSpacing << " [?]" << endl;
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

MDFlexParser::ContainerOption MDFlexParser::getContainerOption() const {
  return containerOption;
}
double MDFlexParser::getCutoff() const { return cutoff; }
MDFlexParser::DataLayoutOption MDFlexParser::getDataLayoutOption() const {
  return dataLayoutOption;
}
MDFlexParser::FunctorOption MDFlexParser::getFunctorOption() const {
  return functorOption;
}

size_t MDFlexParser::getIterations() const {
  return iterations;
}

size_t MDFlexParser::getParticlesPerDim() const {
  return particlesPerDim;
}

double MDFlexParser::getParticlesSpacing() const {
  return particleSpacing;
}