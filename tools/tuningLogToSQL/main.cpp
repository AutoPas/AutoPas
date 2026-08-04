#include "LogToSQLiteWriter.h"

/**
 * This program converts files written by the TuningStrategyLoggerWrapper into a SQLite database file `out.db` for easy
 * data analysis. It already creates some interesting views on the data. See README.md.
 */
int main(int argc, char **argv) {
  if (argc <= 1) {
    std::cerr << "Please provide the data files as arguments" << std::endl;
    exit(1);
  }

  if (not getenv("DISABLE_DEBUG_LOG")) {
    autopas::Logger::get()->set_level(spdlog::level::info);
  }

  LogToSQLiteWriter writer{"out.db"};

  for (int i = 1; i < argc; i++) {
    auto filename = argv[i];
    AutoPasLog(INFO, "Writing file {}: {}", i, filename);

    writer.write(filename);
  }
}
