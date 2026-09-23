#include <algorithm>
#include <iostream>
#include <filesystem>

#include "FilenamesHelper.h"

namespace filenames_helper {
void PrintFilenames(std::vector<std::string> filenames) {
  for (const auto& filename : filenames) {
    std::cout << filename << "\t";
  }
  std::cout << std::endl;
}

std::vector<std::string> GetFilenames(std::string folder_path) {
  std::vector<std::string> filenames;

  if (std::filesystem::is_directory(folder_path)) {
    for (const auto& entry : std::filesystem::directory_iterator(folder_path)) {
      if (std::filesystem::is_regular_file(entry)) {
        std::string filename = entry.path().filename().string();
        if (filename.find("power") != std::string::npos &&
            filename.find("cap") == std::string::npos) {
          filenames.push_back(filename);
        }
      }
    }
  } else {
    std::cerr << "Not a valid directory: " << folder_path << std::endl;
  }

  return filenames;
}

std::vector<std::string> OrderStringsAlphabetically(
    const std::vector<std::string>& input_strings) {
  std::vector<std::string> ordered_strings =
      input_strings;  // Make a copy to preserve the original

  // Use the sort function to order the strings alphabetically
  std::sort(ordered_strings.begin(), ordered_strings.end());

  return ordered_strings;
}

std::vector<std::string> ExtractStringsContainingSubstring(
    const std::vector<std::string>& string_vector,
    const std::string& substring) {
  std::vector<std::string> extracted_strings;

  for (const std::string& str : string_vector) {
    if (str.find(substring) != std::string::npos) {
      extracted_strings.push_back(str);
    }
  }

  return extracted_strings;
}

std::string PopString(std::vector<std::string>& string_vector,
                      const std::string& string_to_pop) {
  auto it =
      std::find(string_vector.begin(), string_vector.end(), string_to_pop);
  if (it != string_vector.end()) {
    std::string popped_value = *it;
    string_vector.erase(it);
    return popped_value;
  } else {
    return "";  // Return an empty string if string_to_pop is not found
  }
}

void HandleFilenamesWithSubstring(std::vector<std::string>& ordered_filenames,
                                  std::vector<std::string>& filenames,
                                  const std::string& substring) {
  std::vector<std::string> matching_filenames =
      ExtractStringsContainingSubstring(filenames, substring);
  std::vector<std::string> ordered_matching_filenames =
      OrderStringsAlphabetically(matching_filenames);
  for (const auto& matching_filename : ordered_matching_filenames) {
    ordered_filenames.push_back(PopString(filenames, matching_filename));
  }
}

std::vector<std::string> ReorderFilenames(std::vector<std::string> filenames) {
  std::vector<std::string> ordered_filenames;

  // Handle "power" filenames
  std::string tmp_string = PopString(filenames, "power");
  if (!tmp_string.empty()) {
    ordered_filenames.push_back(tmp_string);
  }

  // Define the substrings to look for and corresponding sorting priority
  std::vector<std::string> substrings = {"cpu", "memory", "accel"};

  // Handle filenames with each substring
  for (const auto& substring : substrings) {
    HandleFilenamesWithSubstring(ordered_filenames, filenames, substring);
  }

  // Add any remaining filenames to the ordered list
  for (const auto& filename : filenames) {
    ordered_filenames.push_back(filename);
  }

  return ordered_filenames;
}
}  // end namespace filenames_helper