/**
 * @file RDF.h
 * @author D. Martin
 * @date 07.03.2025
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
class RDF {
 public:
  /**
   * Constructor for the RDF
   * @param autoPasContainer: Smart pointer to the Autopas Container
   * @param radiusMin: Minimum radius between particles for the RDF (usually 0)
   * @param radiusMax: Maximum radius between particles for the RDF
   * @param numBins: Number of bins into which the range between radiusMin and radiusMax is divided
   * @param guardArea: Region starting from the simulation border in which particles are not considered as central
   * primary particles for the RDF
   * @param periodicBoundaries: If periodic boundaries should be used
   */
  RDF(const std::shared_ptr<autopas::AutoPas<ParticleType>> autoPasContainer, double radiusMin, double radiusMax,
      size_t numBins, double guardArea, bool periodicBoundaries);

  /**
   * Destructor.
   */
  ~RDF() = default;

  /**
   * Captures the RDF for the current time step. Ideally this is called for every timestep after reaching the
   * equilibrium.
   */
  void captureRDF();

  /**
   * Computes the final RDF as average of all snap-shots capured with captureRDF.
   */
  void computeFinalRDF();

  /**
   * Returns the final RDF
   * @return returns the RDF as pairs of distance and value
   */
  std::vector<std::pair<double, double>> getFinalRDF();

  /**
   * Writes the final RDF to a CSV file
   * @param outputFolder: Output folder to store the CSV file
   * @param filename: name of the csv file
   */
  void writeToCSV(std::string outputFolder, std::string filename);

 protected:
  /**
   * The the nodes' AutoPas container used for simulation.
   */
  std::shared_ptr<autopas::AutoPas<ParticleType>> _autoPasContainer;

  /**
   * Minimum radius between particles for the RDF (usually 0).
   */
  double _radiusMin = 0.0;
  /**
   * Maximum radius between particles for the RDF.
   */
  double _radiusMax = 3.0;
  /**
   * Number of bins into which the range between radiusMin and radiusMax is divided.
   */
  size_t _numBins = 100;
  /**
   * Region starting from the simulation border in which particles are not considered as central
   * primary particles for the RDF.
   */
  double _guardArea = 2.0;
  /**
   * If periodic boundaries should be used.
   */
  bool _periodicBoundaries = false;
  /**
   * Is set to true after computeFinalRDF is called.
   */
  bool _rdfFinished{false};
  /**
   * Stores the snap shots of captureRDF.
   */
  std::vector<std::vector<double>> _rdfs;
  /**
   * Stores the final RDF.
   */
  std::vector<std::pair<double, double>> _finalRdf;
};