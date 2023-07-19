#include "LogToSQLiteWriter.h"

int main(int argc, char** argv) {
  autopas::Logger::create();

  if(argc <= 1) {
    std::cerr << "Please provide the data files as arguments" << std::endl;
    exit(1);
  }

  if(not getenv("DISABLE_DEBUG_LOG")) {
    autopas::Logger::get()->set_level(spdlog::level::info);
  }

  LogToSQLiteWriter writer{"test.db"};

  for(int i = 1; i < argc; i++) {
    auto filename = argv[i];
    AutoPasLog(INFO, "Writing file {}: {}", i, filename);

    writer.write(filename);
  }

  autopas::Logger::unregister();
}