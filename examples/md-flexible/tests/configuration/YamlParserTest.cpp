/**
 * @file YamlParserTest.cpp
 * @author S. Newcome
 * @date 26.06.2026
 */

#include "YamlParserTest.h"

#include <yaml-cpp/yaml.h>

#include "src/configuration/YamlParser.h"
#include "src/domainDecomposition/LoadBalancerOption.h"
#include "testingHelpers/commonTypedefs.h"

using MDFlexParser::YamlParser::parseSequenceOneElementExpected;

/**
 * A scalar node holding a single value is returned unchanged.
 */
TEST_F(YamlParserTest, scalarReturnsValue) {
  const auto node = YAML::Load("InvertedPressure");
  EXPECT_EQ(parseSequenceOneElementExpected(node, "err"), "InvertedPressure");
}

/**
 * A sequence holding exactly one value returns that value (without sequence decoration).
 */
TEST_F(YamlParserTest, singleElementSequenceReturnsValue) {
  const auto node = YAML::Load("[InvertedPressure]");
  EXPECT_EQ(parseSequenceOneElementExpected(node, "err"), "InvertedPressure");
}

/**
 * A sequence with more than one element violates the single-element expectation and throws.
 */
TEST_F(YamlParserTest, multiElementSequenceThrows) {
  const auto node = YAML::Load("[InvertedPressure, None]");
  EXPECT_THROW(parseSequenceOneElementExpected(node, "err"), std::runtime_error);
}

/**
 * By default "all" is rejected: parsed as an option it would expand into every option, violating the
 * single-element expectation.
 */
TEST_F(YamlParserTest, allThrowsByDefault) {
  EXPECT_THROW(parseSequenceOneElementExpected(YAML::Load("all"), "err"), std::runtime_error);
  EXPECT_THROW(parseSequenceOneElementExpected(YAML::Load("[all]"), "err"), std::runtime_error);
}

/**
 * With allThrowsError disabled, "all" is returned verbatim instead of throwing. This is the path the load
 * balancer relies on, where "all" is meant as the single option "A Load-balancing Library (ALL)".
 */
TEST_F(YamlParserTest, allReturnedWhenNotThrowing) {
  EXPECT_EQ(parseSequenceOneElementExpected(YAML::Load("all"), "err", false), "all");
  EXPECT_EQ(parseSequenceOneElementExpected(YAML::Load("[all]"), "err", false), "all");
}

/**
 * End-to-end check of the load balancer special case: the user's "all" is upper-cased to "ALL" and parsed into
 * LoadBalancerOption::all rather than being rejected or expanded into every option.
 */
TEST_F(YamlParserTest, loadBalancerAllParsedAsALL) {
  auto loadBalancerString = parseSequenceOneElementExpected(YAML::Load("[all]"), "err", false);
  if (loadBalancerString == "all") {
    loadBalancerString = "ALL";
  }
  const auto parsedOptions = LoadBalancerOption::parseOptions(loadBalancerString);
  ASSERT_EQ(parsedOptions.size(), 1);
  EXPECT_EQ(*parsedOptions.begin(), LoadBalancerOption::all);
}

/**
 * Tests parsing CubeGrid with particle-density, particle-spacing, and mutual exclusivity.
 */
TEST_F(YamlParserTest, parseCubeGridSpacingAndDensity) {
  MDFlexConfig config;

  // Valid spacing
  {
    const auto node = YAML::Load(
        "particles-per-dimension: [2, 2, 2]\n"
        "particle-spacing: 1.5\n"
        "bottomLeftCorner: [0, 0, 0]\n"
        "particle-type-id: 0\n"
        "velocity: [0, 0, 0]\n");
    std::vector<std::string> errors;
    const auto obj = MDFlexParser::YamlParser::parseCubeGridObject(config, node, errors);
    EXPECT_TRUE(errors.empty());
    EXPECT_DOUBLE_EQ(obj.getParticleSpacing(), 1.5);
    EXPECT_DOUBLE_EQ(obj.getParticleDensity(), 1.0 / (1.5 * 1.5 * 1.5));
  }

  // Valid density
  {
    const auto node = YAML::Load(
        "particles-per-dimension: [2, 2, 2]\n"
        "particle-density: 8.0\n"
        "bottomLeftCorner: [0, 0, 0]\n"
        "particle-type-id: 0\n"
        "velocity: [0, 0, 0]\n");
    std::vector<std::string> errors;
    const auto obj = MDFlexParser::YamlParser::parseCubeGridObject(config, node, errors);
    EXPECT_TRUE(errors.empty());
    EXPECT_DOUBLE_EQ(obj.getParticleSpacing(), 0.5);
    EXPECT_DOUBLE_EQ(obj.getParticleDensity(), 8.0);
  }

  // Both spacing and density specified -> error
  {
    const auto node = YAML::Load(
        "particles-per-dimension: [2, 2, 2]\n"
        "particle-spacing: 1.0\n"
        "particle-density: 1.0\n"
        "bottomLeftCorner: [0, 0, 0]\n"
        "particle-type-id: 0\n"
        "velocity: [0, 0, 0]\n");
    std::vector<std::string> errors;
    MDFlexParser::YamlParser::parseCubeGridObject(config, node, errors);
    EXPECT_FALSE(errors.empty());
  }

  // Neither spacing nor density specified -> error
  {
    const auto node = YAML::Load(
        "particles-per-dimension: [2, 2, 2]\n"
        "bottomLeftCorner: [0, 0, 0]\n"
        "particle-type-id: 0\n"
        "velocity: [0, 0, 0]\n");
    std::vector<std::string> errors;
    MDFlexParser::YamlParser::parseCubeGridObject(config, node, errors);
    EXPECT_FALSE(errors.empty());
  }
}

/**
 * Tests parsing CubeUniform with numberOfParticles, particle-density, and mutual exclusivity.
 */
TEST_F(YamlParserTest, parseCubeUniformCountAndDensity) {
  MDFlexConfig config;

  // Valid numberOfParticles
  {
    const auto node = YAML::Load(
        "numberOfParticles: 100\n"
        "box-length: [2, 2, 2]\n"
        "bottomLeftCorner: [0, 0, 0]\n"
        "particle-type-id: 0\n"
        "velocity: [0, 0, 0]\n");
    std::vector<std::string> errors;
    const auto obj = MDFlexParser::YamlParser::parseCubeUniformObject(config, node, errors);
    EXPECT_TRUE(errors.empty());
    EXPECT_EQ(obj.getParticlesTotal(), 100);
    EXPECT_DOUBLE_EQ(obj.getParticleDensity(), 100.0 / 8.0);
  }

  // Valid density
  {
    const auto node = YAML::Load(
        "particle-density: 2.5\n"
        "box-length: [2, 2, 2]\n"
        "bottomLeftCorner: [0, 0, 0]\n"
        "particle-type-id: 0\n"
        "velocity: [0, 0, 0]\n");
    std::vector<std::string> errors;
    const auto obj = MDFlexParser::YamlParser::parseCubeUniformObject(config, node, errors);
    EXPECT_TRUE(errors.empty());
    EXPECT_EQ(obj.getParticlesTotal(), 20);  // 2.5 * 8 = 20
    EXPECT_DOUBLE_EQ(obj.getParticleDensity(), 2.5);
  }

  // Both count and density specified -> error
  {
    const auto node = YAML::Load(
        "numberOfParticles: 100\n"
        "particle-density: 2.5\n"
        "box-length: [2, 2, 2]\n"
        "bottomLeftCorner: [0, 0, 0]\n"
        "particle-type-id: 0\n"
        "velocity: [0, 0, 0]\n");
    std::vector<std::string> errors;
    MDFlexParser::YamlParser::parseCubeUniformObject(config, node, errors);
    EXPECT_FALSE(errors.empty());
  }

  // Neither count nor density specified -> error
  {
    const auto node = YAML::Load(
        "box-length: [2, 2, 2]\n"
        "bottomLeftCorner: [0, 0, 0]\n"
        "particle-type-id: 0\n"
        "velocity: [0, 0, 0]\n");
    std::vector<std::string> errors;
    MDFlexParser::YamlParser::parseCubeUniformObject(config, node, errors);
    EXPECT_FALSE(errors.empty());
  }
}

/**
 * Tests parsing CubeClosestPacked with structure, particle-density, particle-spacing, and mutual exclusivity.
 */
TEST_F(YamlParserTest, parseCubeClosestPackedStructureAndDensity) {
  MDFlexConfig config;

  // Default structure is FCC, with spacing
  {
    const auto node = YAML::Load(
        "box-length: [4, 4, 4]\n"
        "particle-spacing: 1.0\n"
        "bottomLeftCorner: [0, 0, 0]\n"
        "particle-type-id: 0\n"
        "velocity: [0, 0, 0]\n");
    std::vector<std::string> errors;
    const auto obj = MDFlexParser::YamlParser::parseCubeClosestPacked(config, node, errors);
    EXPECT_TRUE(errors.empty());
    EXPECT_EQ(obj.getLatticeStructure(), CubeClosestPacked::LatticeStructure::FCC);
    EXPECT_DOUBLE_EQ(obj.getParticleSpacing(), 1.0);
  }

  // Explicit HCP structure
  {
    const auto node = YAML::Load(
        "structure: hcp\n"
        "box-length: [4, 4, 4]\n"
        "particle-spacing: 1.0\n"
        "bottomLeftCorner: [0, 0, 0]\n"
        "particle-type-id: 0\n"
        "velocity: [0, 0, 0]\n");
    std::vector<std::string> errors;
    const auto obj = MDFlexParser::YamlParser::parseCubeClosestPacked(config, node, errors);
    EXPECT_TRUE(errors.empty());
    EXPECT_EQ(obj.getLatticeStructure(), CubeClosestPacked::LatticeStructure::HCP);
  }

  // Invalid structure -> error
  {
    const auto node = YAML::Load(
        "structure: invalid_struct\n"
        "box-length: [4, 4, 4]\n"
        "particle-spacing: 1.0\n"
        "bottomLeftCorner: [0, 0, 0]\n"
        "particle-type-id: 0\n"
        "velocity: [0, 0, 0]\n");
    std::vector<std::string> errors;
    MDFlexParser::YamlParser::parseCubeClosestPacked(config, node, errors);
    EXPECT_FALSE(errors.empty());
  }

  // Valid density
  {
    const auto node = YAML::Load(
        "structure: fcc\n"
        "box-length: [4, 4, 4]\n"
        "particle-density: 0.984375\n"
        "bottomLeftCorner: [0, 0, 0]\n"
        "particle-type-id: 0\n"
        "velocity: [0, 0, 0]\n");
    std::vector<std::string> errors;
    const auto obj = MDFlexParser::YamlParser::parseCubeClosestPacked(config, node, errors);
    EXPECT_TRUE(errors.empty());
    EXPECT_EQ(obj.getLatticeStructure(), CubeClosestPacked::LatticeStructure::FCC);
    EXPECT_DOUBLE_EQ(obj.getParticleDensity(), 0.984375);
    EXPECT_DOUBLE_EQ(obj.getParticleSpacing(), std::cbrt(std::sqrt(2.0) / 0.984375));
  }

  // Both spacing and density -> error
  {
    const auto node = YAML::Load(
        "box-length: [4, 4, 4]\n"
        "particle-spacing: 1.0\n"
        "particle-density: 1.0\n"
        "bottomLeftCorner: [0, 0, 0]\n"
        "particle-type-id: 0\n"
        "velocity: [0, 0, 0]\n");
    std::vector<std::string> errors;
    MDFlexParser::YamlParser::parseCubeClosestPacked(config, node, errors);
    EXPECT_FALSE(errors.empty());
  }

  // Neither spacing nor density -> error
  {
    const auto node = YAML::Load(
        "box-length: [4, 4, 4]\n"
        "bottomLeftCorner: [0, 0, 0]\n"
        "particle-type-id: 0\n"
        "velocity: [0, 0, 0]\n");
    std::vector<std::string> errors;
    MDFlexParser::YamlParser::parseCubeClosestPacked(config, node, errors);
    EXPECT_FALSE(errors.empty());
  }

  // Centered: false
  {
    const auto node = YAML::Load(
        "structure: fcc\n"
        "box-length: [4, 4, 4]\n"
        "particle-spacing: 1.0\n"
        "bottomLeftCorner: [1, 2, 3]\n"
        "centered: false\n"
        "particle-type-id: 0\n"
        "velocity: [0, 0, 0]\n");
    std::vector<std::string> errors;
    const auto obj = MDFlexParser::YamlParser::parseCubeClosestPacked(config, node, errors);
    EXPECT_TRUE(errors.empty());
    std::vector<ParticleType> particles;
    obj.generate(particles);
    ASSERT_GT(particles.size(), 0);
    EXPECT_NEAR(particles[0].getR()[0], 1.0, 1e-10);
    EXPECT_NEAR(particles[0].getR()[1], 2.0, 1e-10);
    EXPECT_NEAR(particles[0].getR()[2], 3.0, 1e-10);
  }

  // Centered: true
  {
    const auto node = YAML::Load(
        "structure: fcc\n"
        "box-length: [4, 4, 4]\n"
        "particle-spacing: 1.0\n"
        "bottomLeftCorner: [1, 2, 3]\n"
        "centered: true\n"
        "particle-type-id: 0\n"
        "velocity: [0, 0, 0]\n");
    std::vector<std::string> errors;
    const auto obj = MDFlexParser::YamlParser::parseCubeClosestPacked(config, node, errors);
    EXPECT_TRUE(errors.empty());
    std::vector<ParticleType> particles;
    obj.generate(particles);
    ASSERT_GT(particles.size(), 0);
    const double a = std::sqrt(2.0) * 1.0;
    EXPECT_NEAR(particles[0].getR()[0], 1.0 + a / 4.0, 1e-10);
    EXPECT_NEAR(particles[0].getR()[1], 2.0 + a / 4.0, 1e-10);
    EXPECT_NEAR(particles[0].getR()[2], 3.0 + a / 4.0, 1e-10);
  }
}

/**
 * Tests parsing CubeGrid with centered and uncentered configurations.
 */
TEST_F(YamlParserTest, parseCubeGridAlignment) {
  MDFlexConfig config;

  // centered: false
  {
    const auto node = YAML::Load(
        "particles-per-dimension: [2, 2, 2]\n"
        "particle-spacing: 1.0\n"
        "bottomLeftCorner: [1, 2, 3]\n"
        "centered: false\n"
        "particle-type-id: 0\n"
        "velocity: [0, 0, 0]\n");
    std::vector<std::string> errors;
    const auto obj = MDFlexParser::YamlParser::parseCubeGridObject(config, node, errors);
    EXPECT_TRUE(errors.empty());
    std::vector<ParticleType> particles;
    obj.generate(particles);
    ASSERT_EQ(particles.size(), 8);
    EXPECT_NEAR(particles[0].getR()[0], 1.0, 1e-10);
    EXPECT_NEAR(particles[0].getR()[1], 2.0, 1e-10);
    EXPECT_NEAR(particles[0].getR()[2], 3.0, 1e-10);
  }

  // centered: true
  {
    const auto node = YAML::Load(
        "particles-per-dimension: [2, 2, 2]\n"
        "particle-spacing: 1.0\n"
        "bottomLeftCorner: [1, 2, 3]\n"
        "centered: true\n"
        "particle-type-id: 0\n"
        "velocity: [0, 0, 0]\n");
    std::vector<std::string> errors;
    const auto obj = MDFlexParser::YamlParser::parseCubeGridObject(config, node, errors);
    EXPECT_TRUE(errors.empty());
    std::vector<ParticleType> particles;
    obj.generate(particles);
    ASSERT_EQ(particles.size(), 8);
    EXPECT_NEAR(particles[0].getR()[0], 1.5, 1e-10);
    EXPECT_NEAR(particles[0].getR()[1], 2.5, 1e-10);
    EXPECT_NEAR(particles[0].getR()[2], 3.5, 1e-10);
  }
}