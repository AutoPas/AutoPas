/**
 * @file CLIParser.h
 * @author F. Gratl
 * @date 18.10.2019
 */

#pragma once

#include <getopt.h>

#include <fstream>
#include <iostream>

#include "MDFlexConfig.h"
#include "ParserExitCodes.h"
#include "autopas/utils/TupleUtils.h"

/**
 * Parser for input from the command line.
 */
namespace MDFlexParser::CLIParser {
/**
 * Checks if a checkpoint or yaml file is specified in the given command line arguments.
 * If any files are found, their paths are saves in the respective fields of the given configuration.
 *
 * @param argc number of command line arguments.
 * @param argv command line argument array.
 * @param config configuration where the input is stored.
 */
void inputFilesPresent(int argc, char **argv, MDFlexConfig &config);

/**
 * Parses the Input for the simulation from the command line.
 * @param argc number of command line arguments.
 * @param argv command line argument array.
 * @param config configuration where the input is stored.
 * @return Indicator of success. See MDFlexParser::exitCodes for possible values.
 */
MDFlexParser::exitCodes parseInput(int argc, char **argv, MDFlexConfig &config);

/**
 * Prints the help message to the given stream.
 * @param ostream Typically std::out.
 * @param relPathOfExecutable  Typically argv[0].
 * @param relevantOptions vector of options to include in the help message.
 */
template <class... T>
void printHelpMessage(std::ostream &ostream, const std::string &relPathOfExecutable,
                      const std::tuple<T...> &relevantOptions) {
  // print header
  ostream << "Usage: " << relPathOfExecutable << "\n";
  ostream << "A simple molecular dynamics simulation program showcasing AutoPas.\n\n";
  ostream << "Non-mandatory options to define a simulation:\n";

  // vector that will hold the lines for options
  std::vector<std::string> outputLines;
  // reserve one vector entry for each tuple element
  outputLines.reserve(std::tuple_size_v<std::tuple<T...>>);

  autopas::utils::TupleUtils::for_each(relevantOptions, [&](auto &elem) {
    std::stringstream ss;
    ss << "    --" << std::setw(MDFlexConfig::valueOffset + 2) << std::left << elem.name;
    ss << elem.description;
    ss << '\n';
    outputLines.push_back(ss.str());
  });

  // sort lines alphabetically
  std::sort(std::begin(outputLines), std::end(outputLines));
  // print all lines
  for_each(std::begin(outputLines), std::end(outputLines), [&](auto l) { ostream << l; });

  // print footer
  ostream << '\n';
  ostream << "md-flexible documentation locally available via: 'make doc_doxygen_md-flexible'\n";
  ostream << "and online at: https://autopas.github.io/doxygen_documentation_md-flexible/git-master/\n";
  ostream << "Report bugs to: https://github.com/AutoPas/AutoPas/issues\n";
  ostream << "Full AutoPas documentation at: https://autopas.github.io/doxygen_documentation/git-master/index.html\n";
}

/**
 * Creates a file "_md-flexible" containing the completions definitions for zsh.
 * The function expects possible options to be in the description of an option and in the form of:
 * 'Possible Values: (val1 val2 ... valN)'
 * @param cliOptions vector of options to include in the file.
 */
template <class... T>
void createZSHCompletionFile(const std::tuple<T...> &cliOptions) {
  constexpr auto filename = "_md-flexible";
  std::ofstream fileStream;
  fileStream.open(filename);

  // some header
  fileStream << "#compdef md-flexible\n"
                "\n"
                "typeset -A opt_args\n"
                "local context state line\n"
                "\n"
                "# doc: \n"
                "# THIS FILE IS AUTO GENERATED BY md-flexible --zsh-completions\n"
                "# http://zsh.sourceforge.net/Doc/Release/Completion-System.html#Completion-System\n"
                "# format : \"--arg[description.]: : (list of options)\"\n"
                "\n"
                "_arguments -s -S \\\n";

  // go through all cli options and create the corresponding completion line
  // also, parse possible values from descriptions
  autopas::utils::TupleUtils::for_each(cliOptions, [&](auto &option) {
    // the keyword to look for value suggestions
    constexpr auto keywordForValues{"Possible Values:"};
    constexpr std::pair<char, char> bracketsForValues = std::make_pair('(', ')');
    constexpr auto keywordForPaths{"Path to"};

    // basic name and description for every option.
    fileStream << "     \"--" << option.name << "[" << option.description << "]";

    // If the option requires arguments and the description suggests some, parse them.
    if (option.requiresArgument) {
      std::string choices;
      if (const auto posKeywordForValues = option.description.find(keywordForValues);
          posKeywordForValues != std::string::npos) {
        // start after (+1) first occurrence of '[' after "Possible Values:"
        const auto startPossibleValues{option.description.find(bracketsForValues.first, posKeywordForValues) + 1};
        // end is the first ']' after the start. The length of the substr is the index difference to start.
        const auto strLengthPossibleValues{option.description.find(bracketsForValues.second, startPossibleValues) -
                                           startPossibleValues};
        choices = "(" + option.description.substr(startPossibleValues, strLengthPossibleValues) + ")";
      } else if (const auto posKeywordForPaths = option.description.find(keywordForPaths);
                 posKeywordForPaths != std::string::npos) {
        choices = "_files";
        // look if a filetype is given
        std::regex rgxFileType{"\\.[a-zA-Z]+"};
        std::smatch stringMatch;
        std::regex_search(option.description, stringMatch, rgxFileType);
        if (not stringMatch.empty()) {
          choices += " -g '*" + stringMatch[0].str() + "'";
        }
      }
      // append completion for options
      fileStream << ": :" << choices;
    }

    // closing " and end the line
    fileStream << "\"\\\n";
  });

  fileStream.close();

  std::cout << "Created file: " << filename << std::endl;
}

}  // namespace MDFlexParser::CLIParser
