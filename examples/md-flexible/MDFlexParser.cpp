#include "MDFlexParser.h"

bool MDFlexParser::parseInput(int argc, char **argv) {
  bool displayHelp = false;
  int option, option_index;
  static struct option long_options[] = {
      {"container", required_argument, nullptr, 'c'},
      {"data-layout", required_argument, nullptr, 'd'},
      {"functor", required_argument, nullptr, 'f'},
      {"particles-per-dimension", required_argument, nullptr, 'n'},
  };
  if (argc > 1) {
    string strArg;
    while ((option = getopt_long(argc, argv, "", long_options,
                                 &option_index)) != -1) {
      strArg = optarg;
      transform(strArg.begin(), strArg.end(),
                strArg.begin(), ::tolower);
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
          }
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
          }
          break;
        }
        case 'f': {
          if (strArg.find("lj") != string::npos || strArg.find("lennard-jones") != string::npos) {
            cout << "Using Lennard-Jones (12-6) Functor" << endl;
            functorOption = lj12_6;
          } else {
            cerr << "Unknown functor : " << strArg << endl;
            cerr << "Please use 'Lennard-Jones' or '?'!" << endl;
          }
          break;
        }
        case 'n': {
          try {
            particlesPerDim = (size_t) strtol(optarg, nullptr, 10);
          } catch (const exception &) {
            cerr << "Error parsing number of particles per dimension: " << optarg << endl;
          }
          cout << "Simulating " << particlesPerDim
               << " particles per dimension." << endl;
          break;
        }
        default: {
          // error message handled by getopt
          displayHelp = true;
        }
      }
    }
  } else {
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
MDFlexParser::DataLayoutOption MDFlexParser::getDataLayoutOption() const {
  return dataLayoutOption;
}
size_t MDFlexParser::getParticlesPerDim() const {
  return particlesPerDim;
}

