#include <cstring>
#include <iostream>
#include <iomanip>
#include <memory>
#include <stdexcept>
#include <string>

#include <cxxopts.hpp>

#include <pmt.h>

cxxopts::Options create_commandline_parser(char* argv[]) {
  cxxopts::Options options(argv[0]);

  options.add_options()("n,name", "Name (required)",
                        cxxopts::value<std::string>())(
      "d,device", "Device (optional)",
      cxxopts::value<std::string>()->default_value(""))(
      "command", "Command (optional)",
      cxxopts::value<std::vector<std::string>>()->default_value({}))(
      "h,help", "Print usage");
  options.parse_positional({"command"});

  return options;
}

cxxopts::ParseResult parse_commandline(cxxopts::Options& options, int argc,
                                       char* argv[]) {
  try {
    cxxopts::ParseResult result = options.parse(argc, argv);

    if (result.count("help")) {
      std::cout << options.help() << std::endl;
      exit(EXIT_SUCCESS);
    }

    if (!result.count("name")) {
      throw cxxopts::exceptions::missing_argument("name");
    }

    return result;
  } catch (const cxxopts::exceptions::exception& e) {
    std::cerr << options.help() << std::endl;
    exit(EXIT_FAILURE);
  }
}

void run(pmt::PMT& sensor, const std::vector<std::string>& command) {
  const char* filename = std::getenv(pmt::kDumpFilenameVariable.c_str());
  sensor.StartDump(filename);

  if (command.empty()) {
    while (true) {
      auto state = sensor.Read();
      std::cout << std::fixed << state.timestamp() << ", ";
      for (int i = 0; i < state.NrMeasurements(); i++) {
        std::cout << state.name(i) << ": " << state.watts(i) << " W";
        if (i < (state.NrMeasurements() - 1)) {
          std::cout << ", ";
        }
        std::this_thread::sleep_for(
            std::chrono::milliseconds(sensor.GetMeasurementInterval()));
      }
      std::cout << std::endl;
    }
  } else {
    std::stringstream command_stream;
    for (int i = 1; i < command.size(); i++) {
      if (i > 1) {
        command_stream << " ";
      }
      command_stream << command[i];
    }
    const std::string command_string = command_stream.str();
    auto start = sensor.Read();
    if (system(command_string.c_str()) != 0) {
      perror(command_string.c_str());
    }
    auto end = sensor.Read();
    std::cout << "Runtime: " << pmt::PMT::seconds(start, end) << " s"
              << std::endl;
    std::cout << "Joules: " << pmt::PMT::joules(start, end) << " J"
              << std::endl;
    std::cout << "Watt: " << pmt::PMT::watts(start, end) << " W" << std::endl;
  }
}

int main(int argc, char* argv[]) {
  cxxopts::Options options = create_commandline_parser(argv);
  const cxxopts::ParseResult result = parse_commandline(options, argc, argv);
  const std::string pmt_name = result["name"].as<std::string>();
  const std::string pmt_device = result["device"].as<std::string>();
  std::vector<std::string> command =
      result["command"].as<std::vector<std::string>>();
  if (command.size() == 1 && command[0].empty()) {
    command.resize(0);
  }
  try {
    std::unique_ptr<pmt::PMT> sensor =
        pmt::Create(pmt_name.c_str(), pmt_device.c_str());
    run(*sensor, command);
    return EXIT_SUCCESS;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what();
    return EXIT_FAILURE;
  }
}