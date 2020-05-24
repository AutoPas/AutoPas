/**
 * @file ParserExitCodes.h
 * @author F. Gratl
 * @date 24/05/2020
 */

#pragma once

namespace MDFlexParser {
/**
 * Exit values for the parse functions
 */
enum class exitCodes {
  success,
  parsingError,
  helpFlagFound,
  completionsFlagFound,
};
}