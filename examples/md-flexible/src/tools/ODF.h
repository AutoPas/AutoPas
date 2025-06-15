/**
 * @file ODF.h
 * @author D. Martin
 * @date 12.06.2025
 */
#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "autopas/AutoPasDecl.h"
#include "src/TypeDefinitions.h"

/**
 * Calculates the radial distribution function for a simulation extecuted with md-flexible.
 */
class ODF {
 public:
  /**
   * Constructor for the ODF
   * @param autoPasContainer: Smart pointer to the Autopas Container
   * @param radiusMin: Minimum radius between particles for the ODF (usually 0)
   * @param radiusMax: Maximum radius between particles for the ODF
   * @param numBins: Number of bins into which the range between radiusMin and radiusMax is divided
   * @param guardArea: Region starting from the simulation border in which particles are not considered as central
   * primary particles for the ODF
   * @param periodicBoundaries: If periodic boundaries should be used
   */
  ODF(const std::shared_ptr<autopas::AutoPas<ParticleType>> autoPasContainer, double radiusMin, double radiusMax,
      size_t numBins, double guardArea, bool periodicBoundaries);

  /**
   * Constructor for the ODF that loads data from a CSV file
   * @param filename: Filename to data
   */
  ODF(const std::string &filename);

  /**
   * Destructor.
   */
  ~ODF() = default;

  /**
   * Captures the ODF for the current time step. Ideally this is called for every timestep after reaching the
   * equilibrium.
   */
  void captureODF();

  /**
   * Computes the final ODF as average of all snap-shots capured with captureODF.
   */
  void computeFinalODF();

  /**
   * Returns the final ODF
   * @return returns the ODF as pairs of distance and value
   */
  std::vector<std::pair<double, double>> &getFinalODF();

  /**
   * Writes the final ODF to a CSV file
   * @param outputFolder: Output folder to store the CSV file
   * @param filename: name of the csv file
   */
  void writeToCSV(std::string outputFolder, std::string filename);

  /**
   * Resets the ODF
   */
  void reset();

 protected:
  /**
   * The the nodes' AutoPas container used for simulation.
   */
  std::shared_ptr<autopas::AutoPas<ParticleType>> _autoPasContainer;

  /**
   * Minimum radius between particles for the ODF (usually 0).
   */
  double _radiusMin = 0.0;
  /**
   * Maximum radius between particles for the ODF.
   */
  double _radiusMax = 3.0;
  /**
   * Number of bins into which the range between radiusMin and radiusMax is divided.
   */
  size_t _numBins = 100;
  /**
   * Region starting from the simulation border in which particles are not considered as central
   * primary particles for the ODF.
   */
  double _guardArea = 2.0;
  /**
   * If periodic boundaries should be used.
   */
  bool _periodicBoundaries = false;
  /**
   * Is set to true after computeFinalODF is called.
   */
  bool _odfFinished{false};
  /**
   * Stores the number of captured ODFs
   */
  size_t _numODFs{0};
  /**
   * Stores the final ODF.
   */
  std::vector<std::pair<double, double>> _finalOdf;
};