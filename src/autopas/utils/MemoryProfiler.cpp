/**
 * @file MemoryProfiler.cpp
 * @author F. Gratl
 * @date 9/18/18
 */

#include "MemoryProfiler.h"
#include "Logger.h"

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstring>

size_t autopas::memoryProfiler::currentMemoryUsage() {
  // the stat file should be around ~1400 chars
  static const auto BUFFER_SIZE = 2 * 1024;
  int statusfile = open(statusFileName, O_RDONLY);
  if (statusfile == -1) AutoPasLog(error, "Error opening {}", statusFileName);
  // Advise the kernel of our access pattern.
  posix_fadvise(statusfile, 0, 0, 1);  // FDADVICE_SEQUENTIAL

  char buffer[BUFFER_SIZE + 1];
  // read the file with buffers
  // this loop should only be executed once since the whole file should fit in the buffer.
  while (auto bytes_read = read(statusfile, buffer, BUFFER_SIZE)) {
    if (!bytes_read) break;

    char *needle = buffer;
    // as long as 'V''s are found...
    while (needle != nullptr) {
      // find next start of needle
      needle = (char *)memchr(needle, 'V', (buffer + bytes_read) - needle);
      if (strncmp(needle, "VmRSS:", 6) == 0) {
        // pass whole text after colon. strtol ignores leading whitespaces and stops after number.
        return (size_t)std::strtol(needle + 7, nullptr, 10);
      }
      needle += 6;
    }
  }
  return 0;
}