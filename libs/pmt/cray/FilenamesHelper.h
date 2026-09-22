#ifndef PMT_FILENAMES_HELPER_H
#define PMT_FILENAMES_HELPER_H

#include <string>
#include <vector>

namespace filenames_helper {
// For debug only
void PrintFilenames(std::vector<std::string> filenames);

std::vector<std::string> GetFilenames(std::string folder_path);
std::vector<std::string> OrderStringsAlphabetically(
    const std::vector<std::string>& input_strings);
std::vector<std::string> ExtractStringsContainingSubstring(
    const std::vector<std::string>& string_vector,
    const std::string& substring);
std::string PopString(std::vector<std::string>& string_vector,
                      const std::string& string_to_pop);
void HandleFilenamesWithSubstring(std::vector<std::string>& ordered_filenames,
                                  std::vector<std::string>& filenames,
                                  const std::string& substring);
std::vector<std::string> ReorderFilenames(std::vector<std::string> filenames);

}  // end namespace filenames_helper

#endif
