/**
 * @file OpenMPChunkSizeTest.cpp
 * @author S. J. Newcome
 * @date 13/10/2026
 */

#include <gtest/gtest.h>

#include <set>
#include <sstream>
#include <string>

#include "autopas/utils/OpenMPChunkSize.h"

using namespace autopas;

/**
 * Checks that to_string() and parse() are the inverse of each other for all types of chunk sizes, and that parse()
 * rejects everything which is not a chunk size.
 */
TEST(OpenMPChunkSizeTest, testToStringAndParse) {
  EXPECT_EQ(OpenMPChunkSize(4).to_string(), "4");
  EXPECT_EQ(OpenMPChunkSize::balanced.to_string(), "balanced");
  EXPECT_EQ(OpenMPChunkSize::golden.to_string(), "golden");

  EXPECT_EQ(OpenMPChunkSize::parse("4"), OpenMPChunkSize(4));
  EXPECT_EQ(OpenMPChunkSize::parse("balanced"), OpenMPChunkSize::balanced);
  EXPECT_EQ(OpenMPChunkSize::parse("golden"), OpenMPChunkSize::golden);
  EXPECT_EQ(OpenMPChunkSize::parse("Golden"), OpenMPChunkSize::golden) << "parse() should ignore the case.";

  EXPECT_ANY_THROW(OpenMPChunkSize::parse(""));
  EXPECT_ANY_THROW(OpenMPChunkSize::parse("0")) << "A fixed chunk size has to be at least 1.";
  EXPECT_ANY_THROW(OpenMPChunkSize::parse("-4"));
  EXPECT_ANY_THROW(OpenMPChunkSize::parse("4.5"));
  EXPECT_ANY_THROW(OpenMPChunkSize::parse("4x"));
  EXPECT_ANY_THROW(OpenMPChunkSize::parse("silver"));
}

/**
 * Checks that parseChunkSizes() parses a list of fixed and special chunk sizes with different delimiters, and rejects
 * the whole list if one entry is not a chunk size.
 */
TEST(OpenMPChunkSizeTest, testParseChunkSizes) {
  const std::set<OpenMPChunkSize> expectedChunkSizes{1, 16, OpenMPChunkSize::balanced, OpenMPChunkSize::golden};
  EXPECT_EQ(OpenMPChunkSize::parseChunkSizes("1, 16, balanced, golden"), expectedChunkSizes);
  EXPECT_EQ(OpenMPChunkSize::parseChunkSizes("golden;16 balanced|1"), expectedChunkSizes);

  EXPECT_TRUE(OpenMPChunkSize::parseChunkSizes("").empty());
  EXPECT_ANY_THROW(OpenMPChunkSize::parseChunkSizes("1, silver, golden"));
}

/**
 * Checks the comparison operators: chunk sizes are only equal to themselves, fixed chunk sizes are ordered by size and
 * come before balanced, which comes before golden.
 */
TEST(OpenMPChunkSizeTest, testComparison) {
  EXPECT_EQ(OpenMPChunkSize(4), OpenMPChunkSize(4));
  EXPECT_NE(OpenMPChunkSize(4), OpenMPChunkSize(5));
  EXPECT_NE(OpenMPChunkSize::balanced, OpenMPChunkSize::golden);
  EXPECT_NE(OpenMPChunkSize::balanced, OpenMPChunkSize(0));

  EXPECT_LT(OpenMPChunkSize(4), OpenMPChunkSize(5));
  EXPECT_LT(OpenMPChunkSize(1000), OpenMPChunkSize::balanced);
  EXPECT_LT(OpenMPChunkSize::balanced, OpenMPChunkSize::golden);
}

/**
 * Checks that equal chunk sizes have the same hash, and that the hashes of some different chunk sizes differ.
 */
TEST(OpenMPChunkSizeTest, testHash) {
  const std::hash<OpenMPChunkSize> hash{};
  EXPECT_EQ(hash(OpenMPChunkSize(4)), hash(OpenMPChunkSize(4)));
  EXPECT_EQ(hash(OpenMPChunkSize::golden), hash(OpenMPChunkSize::golden));

  EXPECT_NE(hash(OpenMPChunkSize(4)), hash(OpenMPChunkSize(5)));
  EXPECT_NE(hash(OpenMPChunkSize::balanced), hash(OpenMPChunkSize::golden));
  EXPECT_NE(hash(OpenMPChunkSize::balanced), hash(OpenMPChunkSize(0)));
}

/**
 * Checks that the stream operators write and read chunk sizes, and that reading stops at the end of the chunk size and
 * leaves the rest in the stream.
 */
TEST(OpenMPChunkSizeTest, testStreamOperators) {
  std::stringstream stream;
  stream << OpenMPChunkSize(4) << " , " << OpenMPChunkSize::golden << "}";
  EXPECT_EQ(stream.str(), "4 , golden}");

  OpenMPChunkSize chunkSize;
  stream >> chunkSize;
  EXPECT_EQ(chunkSize, OpenMPChunkSize(4));
  std::string separator;
  stream >> separator;
  EXPECT_EQ(separator, ",");
  stream >> chunkSize;
  EXPECT_EQ(chunkSize, OpenMPChunkSize::golden);
  EXPECT_EQ(stream.get(), '}');
}