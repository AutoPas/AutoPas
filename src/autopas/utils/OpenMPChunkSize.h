/**
 * @file OpenMPChunkSize.h
 * @author S. J. Newcome
 * @date 12/10/2026
 */

#pragma once

#include <algorithm>
#include <cctype>
#include <compare>
#include <cstddef>
#include <functional>
#include <istream>
#include <ostream>
#include <set>
#include <string>

#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/HashCombine.h"
#include "autopas/utils/StringUtils.h"
#include "autopas/utils/logging/Logger.h"

namespace autopas {

/**
 * Class representing the choices for OpenMP's scheduling chunk size: a fixed size (i.e. a simple size_t), or one of the
 * special chunk sizes (balanced and golden).
 *
 * Every object stores a Type and a size. The size is only used if the type is fixed, and is 0 otherwise.
 * - A fixed chunk size OpenMPChunkSize is constructed from a size_t, e.g. `OpenMPChunkSize chunkSize = 4;`.
 * - The special chunk sizes are not constructed. They are the constants OpenMPChunkSize::balanced and
 *   OpenMPChunkSize::golden, which are used like enumerators,
 *   e.g. `OpenMPChunkSize chunkSize = OpenMPChunkSize::golden`.
 */
class OpenMPChunkSize {
 public:
  /**
   * Possible types of chunk sizes.
   */
  enum class Type {
    /**
     * A fixed number of iterations per chunk.
     */
    fixed,
    /**
     * OpenMP's default for static scheduling: the iterations are divided into one chunk per thread, of (approximately)
     * equal size. Only compatible with OpenMPKindOption::omp_static.
     */
    balanced,
    /**
     * LB4OMP's golden (in papers, referred to as expert) chunk size, which depends on the number of iterations and
     * threads of the loop (see OpenMPConfigurator::getGoldenChunkSize()).
     */
    golden,
  };

  /**
   * Constructor. Creates a fixed chunk size of 0, which is not a valid chunk size.
   */
  constexpr OpenMPChunkSize() = default;

  /**
   * Constructor for a fixed chunk size.
   * @param size Number of iterations per chunk.
   */
  // NOLINTNEXTLINE: Not explicit, so that a number can be given wherever a chunk size is expected.
  constexpr OpenMPChunkSize(size_t size) : _size(size) {}

  /**
   * The balanced chunk size.
   */
  static const OpenMPChunkSize balanced;

  /**
   * The golden chunk size.
   */
  static const OpenMPChunkSize golden;

  /**
   * Getter for the type of the chunk size.
   * @return
   */
  [[nodiscard]] constexpr Type getType() const { return _type; }

  /**
   * Getter for the number of iterations per chunk. Only meaningful for fixed chunk sizes, warns otherwise.
   * @return The size for fixed chunk sizes, otherwise 0.
   */
  [[nodiscard]] size_t getSize() const {
    if (_type != Type::fixed) {
      AutoPasLog(WARN, "OpenMPChunkSize::getSize(): The chunk size {} is not a fixed chunk size, so it has no size.",
                 to_string());
    }
    return _size;
  }

  /**
   * Comparison operator. Fixed chunk sizes are ordered by size and come before balanced, which comes before golden.
   * @return
   */
  constexpr auto operator<=>(const OpenMPChunkSize &) const = default;

  /**
   * Provides a string representation of the chunk size.
   * @return The size for fixed chunk sizes, otherwise the name of the special chunk size.
   */
  [[nodiscard]] std::string to_string() const {
    switch (_type) {
      case Type::fixed:
        return std::to_string(_size);
      case Type::balanced:
        return "balanced";
      case Type::golden:
        return "golden";
    }
    utils::ExceptionHandler::exception("OpenMPChunkSize::to_string(): Unknown type {}.", static_cast<int>(_type));
    return "";
  }

  /**
   * Converts a string to a chunk size. Inverse of to_string(), ignoring case.
   * @param chunkSizeString A positive integer, "balanced" or "golden".
   * @return The chunk size.
   */
  static OpenMPChunkSize parse(const std::string &chunkSizeString) {
    std::string lowerCaseString = chunkSizeString;
    std::transform(lowerCaseString.begin(), lowerCaseString.end(), lowerCaseString.begin(), ::tolower);
    if (lowerCaseString == "balanced") {
      return balanced;
    }
    if (lowerCaseString == "golden") {
      return golden;
    }

    // Check that a non-balanced/golden chunk size is actually a positive integer value
    if (lowerCaseString.empty() or not std::all_of(lowerCaseString.begin(), lowerCaseString.end(), ::isdigit)) {
      utils::ExceptionHandler::exception(
          "OpenMPChunkSize::parse(): \"{}\" is not a chunk size. Possible values: a positive integer, balanced, "
          "golden.",
          chunkSizeString);
    }
    const auto size = utils::StringUtils::parseNumber<size_t>(lowerCaseString);
    if (size == 0) {
      utils::ExceptionHandler::exception(
          "OpenMPChunkSize::parse(): A fixed chunk size must be at least 1. For a default chunk size for static, use "
          "'balanced'.");
    }
    return {size};
  }

  /**
   * Converts a string of chunk sizes to a set of chunk sizes.
   * @param chunkSizesString Chunk sizes as parse() takes them, separated by autopas::utils::StringUtils::delimiters.
   * @return Set of chunk sizes. Empty if the string contains no chunk sizes.
   */
  static std::set<OpenMPChunkSize> parseChunkSizes(const std::string &chunkSizesString) {
    std::set<OpenMPChunkSize> chunkSizes;
    for (const auto &token : utils::StringUtils::tokenize(chunkSizesString, utils::StringUtils::delimiters)) {
      chunkSizes.insert(parse(token));
    }
    return chunkSizes;
  }

  /**
   * Stream output operator.
   * @param os
   * @param chunkSize
   * @return
   */
  friend std::ostream &operator<<(std::ostream &os, const OpenMPChunkSize &chunkSize) {
    os << chunkSize.to_string();
    return os;
  }

  /**
   * Stream extraction operator.
   * @param in Stream from which will be parsed.
   * @param chunkSize out-parameter that will be filled with the parsed chunk size.
   * @return Processed stream.
   */
  friend std::istream &operator>>(std::istream &in, OpenMPChunkSize &chunkSize) {
    in >> std::ws;
    // A stream does not say where the chunk size ends, e.g. in "4 , VectorizationPattern: ...". So take everything up
    // to the first character which cannot be part of a chunk size, and leave the rest in the stream.
    std::string str;
    while (std::isalnum(in.peek())) {
      str.push_back(static_cast<char>(in.get()));
    }
    chunkSize = parse(str);
    return in;
  }

 private:
  /**
   * Constructor for the special chunk sizes. Only used for the definitions of balanced and golden below the class.
   * Fixed chunk sizes are constructed from their size instead.
   * @param type Type::balanced or Type::golden.
   */
  constexpr explicit OpenMPChunkSize(Type type) : _type(type) {
    if (type == Type::fixed) {
      utils::ExceptionHandler::exception("OpenMPChunkSize: A fixed chunk size has to be constructed from its size.");
    }
  }

  /**
   * The type of the chunk size.
   */
  Type _type{Type::fixed};

  /**
   * Number of iterations per chunk.
   */
  size_t _size{0};
};

// Expose balanced and golden similarly to if they were enumerators. Under the hood, they are objects of OpenMPChunkSize
// with the type balanced and golden respectively, created using the private constructor above. They are defined here
// and not inside the class, as a class cannot contain objects of itself before its definition is complete.
inline constexpr OpenMPChunkSize OpenMPChunkSize::balanced{OpenMPChunkSize::Type::balanced};
inline constexpr OpenMPChunkSize OpenMPChunkSize::golden{OpenMPChunkSize::Type::golden};

}  // namespace autopas

/**
 * Hash for OpenMPChunkSize objects, to be used in e.g. unordered maps or in the hash of a Configuration.
 */
template <>
struct std::hash<autopas::OpenMPChunkSize> {
  /**
   * Hash Function operator
   * @param chunkSize
   * @return
   */
  std::size_t operator()(const autopas::OpenMPChunkSize &chunkSize) const {
    // Only fixed chunk sizes have a size.
    const size_t size = chunkSize.getType() == autopas::OpenMPChunkSize::Type::fixed ? chunkSize.getSize() : 0;
    return autopas::utils::hashCombine(chunkSize.getType(), size);
  }
};