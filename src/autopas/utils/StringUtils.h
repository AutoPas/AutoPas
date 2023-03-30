/**
 * @file StringUtils.h
 * @author F. Gratl
 * @date 1/14/19
 */

#pragma once

#include <cmath>
#include <regex>
#include <set>
#include <string>
#include <vector>

#include "autopas/utils/NumberInterval.h"
#include "autopas/utils/NumberSet.h"
#include "autopas/utils/NumberSetFinite.h"

namespace autopas::utils::StringUtils {

// anonymous namespace for namespace-private helper functions
namespace {
/**
 * Calculates a similarity score of s1 and s2 based on the Needleman-Wunsch string alignment algorithm.
 *
 * @param s1
 * @param s2
 * @return Score in the lower right corner of the score matrix.
 */
inline int needlemanWunschScore(std::string s1, std::string s2) {
  // these scores correspond to the number of edits needed to match s1 to s2
  constexpr int scoreMatch = 1;
  constexpr int scoreMismatch = -1;
  constexpr int scoreGap = -1;

  // |s1|+1 x |s2|+1 Matrix
  std::vector<std::vector<int>> scoreMatrix(s1.length() + 1, std::vector<int>(s2.length() + 1, 0));

  // initialize top and right border with cumulative gap penalties
  for (size_t i = 0; i < scoreMatrix.size(); ++i) {
    scoreMatrix[i][0] = i * scoreGap;
  }
  for (size_t j = 0; j < scoreMatrix[0].size(); ++j) {
    scoreMatrix[0][j] = j * scoreGap;
  }

  // fill rest of matrix
  for (size_t i = 1; i < scoreMatrix.size(); ++i) {
    for (size_t j = 1; j < scoreMatrix[0].size(); ++j) {
      auto matchValue = s1[i - 1] == s2[j - 1] ? scoreMatch : scoreMismatch;
      auto scoreDiagonal = scoreMatrix[i - 1][j - 1] + matchValue;
      auto scoreLeft = scoreMatrix[i - 1][j] + scoreGap;
      auto scoreTop = scoreMatrix[i][j - 1] + scoreGap;

      std::array<decltype(scoreDiagonal), 3> scores = {scoreDiagonal, scoreLeft, scoreTop};
      auto scoreMax = std::max_element(scores.begin(), scores.end());

      scoreMatrix[i][j] = *scoreMax;
    }
  }

  // omit backtracking since we are not interested in the alignment but only in
  // the score lower right corner contains similarity score
  return scoreMatrix[scoreMatrix.size() - 1][scoreMatrix[scoreMatrix.size() - 1].size() - 1];
}
}  // namespace
/**
 * Finds best match of needle in haystack.
 *
 * Needle is compared to every option in haystack and the Needleman-Wunsch score calculated.
 * If the result is ambiguous an exception is thrown.
 *
 * @param haystack Vector of string to match to.
 * @param needle
 * @return Best matching string.
 */
inline std::string matchStrings(const std::vector<std::string> &haystack, std::string needle) {
  std::transform(needle.begin(), needle.end(), needle.begin(), ::tolower);
  auto bestDistance = std::numeric_limits<int>::min();
  std::vector<std::string> matchedStrings;
  for (auto &s : haystack) {
    auto distance = needlemanWunschScore(needle, s);
    // if we find a better match throw out current matches
    if (distance > bestDistance) {
      matchedStrings.clear();
      bestDistance = distance;
    }
    // save every match that is at least as good as the current one
    if (distance >= bestDistance) {
      matchedStrings.push_back(s);
    }
  }
  if (matchedStrings.size() > 1) {
    utils::ExceptionHandler::exception("Given String ({}) is ambiguous! Which option do you mean: {}", needle,
                                       [](auto arr) -> std::string {
                                         std::ostringstream ss;
                                         for (auto &a : arr) {
                                           ss << a << ", ";
                                         }
                                         // deletes last comma
                                         ss << "\b\b";
                                         return ss.str();
                                       }(matchedStrings));
  }
  return matchedStrings[0];
}

/**
 * All accepted delimiters to split input strings.
 */
constexpr char delimiters[] = " ,;|/";
/**
 * Regex for all delimiters to split input strings.
 */
constexpr char delimitersRgx[] = "[\\s,;|/]";
/**
 * Regex for all but delimiters to split input strings as regex.
 */
constexpr char delimitersRgxInv[] = "[^\\s,;|/]";

/**
 *  Regex for a double e.g. 1 | 1.2 | 1.2e-3
 */
static const std::string regexDoubleStr{
    "[0-9]+"  // at least one int
    "\\.?"    // maybe a dot
    "[0-9]*"  // maybe more integers after the dot
    "(?:"     // start of non-capturing group for exp
    "e"       // exponent
    "-?"      // optional minus
    "[0-9]+"  // at least one int
    ")?"      // end of group, group is optional
};


/**
 *  Regex for a int e.g. 1 | 2 | 3
 */
static const std::string regexIntStr{
    "[0-9]+"  // at least one int
};

/**
 * Splits a string by multiple delimiters.
 * @param searchString
 * @param delimiters
 * @return Set of substrings.
 */
inline std::vector<std::string> tokenize(const std::string &searchString, const std::string &delimiters) {
  std::vector<std::string> wordVector;

  std::size_t prev = 0, pos;
  while ((pos = searchString.find_first_of(delimiters, prev)) != std::string::npos) {
    if (pos > prev) wordVector.push_back(searchString.substr(prev, pos - prev));
    prev = pos + 1;
  }
  if (prev < searchString.length()) wordVector.push_back(searchString.substr(prev, std::string::npos));

  return wordVector;
}

/**
 * Converts a string to std::array<double,3>.
 * Allowed delimiters can be found in autopas::utils::StringUtils::delimiters.
 *
 * String format: 3 doubles(or ints) separated by delimiters (examples: 10.,10.,10.)
 *
 * @param string String to parse.
 * @return
 */
inline std::array<double, 3> parseArrayD3(const std::string &string) {
  std::array<double, 3> parsedArray{};
  auto strings = tokenize(string, delimiters);
  if (strings.size() > 3) {
    autopas::utils::ExceptionHandler::exception("parseArrayD3(): found {} instead of 3 array fields.", strings.size());
  }
  for (int i = 0; i < 3; i++) {
    try {
      parsedArray[i] = std::stod(strings[i]);
    } catch (const std::exception &e) {
      autopas::utils::ExceptionHandler::exception("parseArrayD3(): could not convert {} to a double: \n{}", strings[i],
                                                  e.what());
    }
  }
  return parsedArray;
}

/**
 * Converts a string to bool
 *
 * String format: on || off || enabled || disabled || true || false
 * @param booleanOption
 * @return
 */
inline bool parseBoolOption(const std::string &booleanOption) {
  if (booleanOption == "on" or booleanOption == "true" or booleanOption == "enabled") {
    return true;
  } else if (booleanOption == "off" or booleanOption == "false" or booleanOption == "disabled") {
    return false;
  } else {
    autopas::utils::ExceptionHandler::exception("Unknown boolean Option: {}", booleanOption);
  }
  // should not be reached
  return false;
}

/**
 * Converts a string to a set of doubles.
 * @param doubleString String containing doubles.
 * @return Set of doubles. If no valid double was found the empty set is returned.
 */
inline std::set<double> parseDoubles(const std::string &doubleString) {
  std::set<double> doubles;

  std::regex regexDouble(regexDoubleStr);

  // use regex iter to find all doubles in the string.
  for (auto number = std::sregex_iterator(doubleString.begin(), doubleString.end(), regexDouble);
       number != std::sregex_iterator(); ++number) {
    try {
      double value = stod(number->str());
      doubles.insert(value);
    } catch (const std::exception &) {
      autopas::utils::ExceptionHandler::exception("Failed to parse a double from: {}", number->str());
    }
  }

  return doubles;
}

/**
 * Converts a string to a set of ints.
 * @param intString String containing ints.
 * @return Set of ints. If no valid int was found the empty set is returned.
 */
inline std::set<int> parseInts(const std::string &intString) {
  std::set<int> ints;

  std::regex regexInt(regexIntStr);

  // use regex iter to find all ints in the string.
  for (auto number = std::sregex_iterator(intString.begin(), intString.end(), regexInt);
       number != std::sregex_iterator(); ++number) {
    try {
      int value = stoi(number->str());
      ints.insert(value);
    } catch (const std::exception &) {
      autopas::utils::ExceptionHandler::exception("Failed to parse a int from: {}", number->str());
    }
  }

  return ints;
}

/**
 * Converts a string to a NumberSet<double>.
 *
 * @note Formats:
 * NumberSetFinite [x,y,z]
 * NumberInterval x-y
 *
 * @param setString String containing the set.
 * @return NumberSet<double>. If no valid double was found the empty set is returned.
 */
inline std::unique_ptr<autopas::NumberSet<double>> parseNumberSetDoubles(const std::string &setString) {
  // try to match an interval x-y
  std::regex regexInterval("("                 // start of 1. capture
                           + regexDoubleStr +  // a double
                           ")"                 // end of 1. capture
                           "\\s*"              // maybe whitespaces
                           "-"                 // a dash
                           "\\s*"              // maybe more whitespaces
                           "("                 // start of 2. capture
                           + regexDoubleStr +  // a double
                           ")"                 // end of 2. capture
  );
  std::smatch matches ;
  if (std::regex_match(setString, matches, regexInterval)) {
    try {
      // matchers has whole string as str(0) so start at 1
      double min = stod(matches.str(1));
      double max = stod(matches.str(2));
      return std::make_unique<autopas::NumberInterval<double>>(min, max);
    } catch (const std::exception &) {
      // try parseDoubles instead
    }
  }

  std::set<double> values =   autopas::utils::StringUtils::parseDoubles(setString);
  return std::make_unique<autopas::NumberSetFinite<double>>(values);
}

/**
 * Converts a string to a NumberSet<int>.
 *
 * @note Formats:
 * NumberSetFinite [x,y,z]
 * NumberInterval x-y
 *
 * @param setString String containing the set.
 * @return NumberSet<int>. If no valid int was found the empty set is returned.
 */
inline std::unique_ptr<autopas::NumberSet<int>> parseNumberSetInts(const std::string &setString) {
  // try to match an interval x-y
  std::regex regexInterval("("                 // start of 1. capture
                           + regexIntStr +  // a double
                           ")"                 // end of 1. capture
                           "\\s*"              // maybe whitespaces
                           "-"                 // a dash
                           "\\s*"              // maybe more whitespaces
                           "("                 // start of 2. capture
                           + regexIntStr +  // a double
                           ")"                 // end of 2. capture
  );
  std::smatch matches;
  if (std::regex_match(setString, matches, regexInterval)) {
    try {
      std::set<int> numbers =  std::set<int>({});
      for (const auto& match : matches) {
        numbers.insert(std::stoi(match.str()));
      }
      return std::make_unique<autopas::NumberSetFinite<int>>(std::set<int>({numbers}));
    } catch (const std::exception &) {
      // try parseDoubles instead
    }
  }

  std::set<int> values = autopas::utils::StringUtils::parseInts(setString);
  return std::make_unique<autopas::NumberSetFinite<int>>(values);
}

}  // namespace autopas::utils::StringUtils