#pragma once

#include "autopas/selectors/Configuration.h"
#include "autopas/selectors/tuningStrategy/LiveInfo.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <tuple>

#include "sqlite3.h"

namespace autopas::tuningLogEntry {
template <class T>
void toStringHelper(std::ostream &in, const T &val) {
  in << val << ' ';
}

template <class... Payload>
auto toString(const Payload &... payload) {
  std::stringstream stream;
  (toStringHelper(stream, payload), ...);
  return stream.str();
}

template <class... Payload>
std::tuple<Payload...> fromString(std::stringstream &stream) {
  std::tuple<Payload...> tuple{};
  ((stream >> std::get<Payload>(tuple)), ...);
  return tuple;
}

std::string writeEvidence(long time, size_t iteration, const Configuration &config) {
  return toString(std::string{"evidence"}, time, iteration, config);
}

std::tuple<long, size_t, Configuration> readEvidence(std::stringstream &str) {
  return fromString<long, size_t, Configuration>(str);
}

std::string writeTune(bool currentInvalid) { return toString(std::string{"tune"}, currentInvalid); }

bool readTune(std::stringstream &str) { return std::get<0>(fromString<bool>(str)); }

std::string writeReset(size_t iteration) { return toString(std::string{"reset"}, iteration); }

size_t readReset(std::stringstream &str) { return std::get<0>(fromString<size_t>(str)); }

std::string writeLiveInfo(const LiveInfo &liveInfo) { return toString(std::string{"liveInfo"}, liveInfo); }

LiveInfo readLiveInfo(std::stringstream &str) { return std::get<0>(fromString<LiveInfo>(str)); }
};  // namespace autopas::tuningLogEntry

static int callback(void*, int argc, char** argv, char** azColName) {
  for(int i = 0; i < argc; ++i) {
    std::cout << azColName[i] << " = " << (argv[i] ? argv[i] : "NULL") << " ";
  }
  std::cout << std::endl;
  return 0;
}

static std::string escapeSingleQuotesForSQL(const std::string& str) {
  std::string result;
  result.reserve(str.size() + 6);

  for(char c : str) {
    result.push_back(c);
    if(c == '\'') {
      result.push_back('\'');
    }
  }

  return result;
}

class LogToSQLiteWriter {
 public:
  explicit LogToSQLiteWriter(const std::string& databaseName) : _db(nullptr) {
    auto failed = sqlite3_open(databaseName.c_str(), &_db);

    if(failed) {
      AutoPasLog(error, "Can't open database {}: {}", databaseName, sqlite3_errmsg(_db));
      throw std::invalid_argument{databaseName};
    }

    createSchema();
  }

  void write(const char* filename) {
    std::ifstream in{filename};

    if (not in.is_open()) {
      AutoPasLog(error, "Could not open file {}", filename);
      AutoPasLog(error, "Exiting!");
      exit(-1);
    }

    auto escapedFilename = escapeSingleQuotesForSQL(filename);
    std::stringstream insertScenarios;
    std::stringstream insertMeasurements;
    insertScenarios << "INSERT INTO ScenarioRaw VALUES ";
    insertMeasurements << "INSERT INTO MeasurementRaw VALUES ";
    char sepScenarios = ' ';
    char sepMeasurements = ' ';
    while (not in.eof()) {
      std::string line;
      std::getline(in, line, '\n');

      std::stringstream stream{line};
      std::string type;
      std::getline(stream, type, ' ');

      if (type == "evidence") {
        const auto &[time, iteration, config] = autopas::tuningLogEntry::readEvidence(stream);
        insertMeasurements << sepMeasurements << "(\'" << escapedFilename << "\',\'"
            << config.container << "\'," << config.cellSizeFactor << ",\'" << config.traversal << "\',\'"
            << config.loadEstimator << "\',\'" << config.dataLayout << "\',\'" << config.newton3 << "\',"
            << iteration << "," << time << ")";
        sepMeasurements = ',';
      } else if (type == "tune") {
        // Do nothing in former tune
      } else if (type == "liveInfo") {
        const auto &liveInfo = autopas::tuningLogEntry::readLiveInfo(stream);
        const auto& d = liveInfo.get();
        auto toStr = [](const auto& variant) {
          return std::visit([](const auto &val) { return std::to_string(val); }, variant);
        };
        insertScenarios << sepScenarios << "( \'" << escapedFilename << "\',"
            << toStr(d.at("avgParticlesPerCell")) << "," << toStr(d.at("cutoff")) << ","
            << toStr(d.at("domainSizeX")) << "," << toStr(d.at("domainSizeY")) << ","
            << toStr(d.at("domainSizeZ")) << ","
                        << toStr(d.at("maxParticlesPerCell")) << ","
            << toStr(d.at("minParticlesPerCell")) << "," << toStr(d.at("numCells")) << ","
            << toStr(d.at("numEmptyCells")) << "," << toStr(d.at("numHaloParticles")) << ","
            << toStr(d.at("numParticles")) << "," << toStr(d.at("particleSize")) << ","
            << toStr(d.at("particleSizeNeededByFunctor")) << ","
            << toStr(d.at("particlesPerBlurredCellStdDev")) << ","
            << toStr(d.at("particlesPerCellStdDev")) << ","
            << toStr(d.at("skin")) << "," << toStr(d.at("threadCount")) << ","
            << toStr(d.at("rebuildFrequency")) << "," << toStr(d.at("estimatedNumNeighborInteractions"))
            << ")";
        sepScenarios = ',';
      } else if (type == "reset" || in.eof()) {
        // Do nothing on reset
      }
    }

    std::cout << insertScenarios.str() << std::endl;
    std::cout << insertMeasurements.str() << std::endl;
    sendQuery(insertScenarios.str().c_str());
    sendQuery(insertMeasurements.str().c_str());
  }

 ~LogToSQLiteWriter() {
   sqlite3_close(_db);
 }

private:
 void createSchema() {
   sendQuery(schemaSQL);
 }

 void sendQuery(const char* query) {
   char* errMsg = nullptr;
   auto returnCode = sqlite3_exec(_db, query, callback, nullptr, &errMsg);

   if(returnCode != SQLITE_OK) {
     AutoPasLog(error, "SQL error: {}", errMsg);
     sqlite3_free(errMsg);
   }
 }

 private:
  sqlite3* _db;

  static constexpr const char* schemaSQL = R"SQL(
CREATE TABLE ScenarioRaw (
  filename VARCHAR(1024) PRIMARY KEY,
  avgParticlesPerCell FLOAT,
  cutoff FLOAT,
  domainSizeX FLOAT,
  domainSizeY FLOAT,
  domainSizeZ FLOAT,
  maxParticlesPerCell INTEGER,
  minParticlesPerCell INTEGER,
  numCells INTEGER,
  numEmptyCells INTEGER,
  numHaloParticles INTEGER,
  numParticles INTEGER,
  particleSize INTEGER,
  particleSizeNeededByFunctor INTEGER,
  particlesPerBlurredCellStdDev FLOAT,
  particlesPerCellStdDev FLOAT,
  skin FLOAT,
  threadCount INTEGER,
  rebuildFrequency INTEGER,
  estimatedNumNeighborInteractions INTEGER
);
CREATE TABLE MeasurementRaw (
  scenario VARCHAR(1024),
  container VARCHAR(64),
  cellSizeFactor FLOAT,
  traversal VARCHAR(64),
  loadEstimator VARCHAR(64),
  dataLayout VARCHAR(64),
  newton3 VARCHAR(64),
  iteration INTEGER,
  nsRuntime BIGINT,
  FOREIGN KEY (scenario) REFERENCES Scenario(filename),
  PRIMARY KEY (
    scenario,
    container,
    cellSizeFactor,
    traversal,
    loadEstimator,
    dataLayout,
    newton3
  )
);

CREATE VIEW Scenario AS
SELECT
  *
FROM
  ScenarioRaw
;

CREATE VIEW Measurement AS
SELECT
  *
FROM
  MeasurementRaw
;

CREATE VIEW numDifferentConfigs AS
SELECT
  COUNT(*)
FROM
  (SELECT
     1
   FROM
     Measurement
   GROUP BY container, traversal, dataLayout, newton3, loadEstimator
  )
;

CREATE VIEW bestConfigurations AS
SELECT
  container,
  dataLayout,
  newton3,
  loadEstimator,
  traversal,
  MIN(nsRuntime) as nsMinRuntime,
  scenario
FROM
  Measurement
group by
  scenario
ORDER BY
  container,
  dataLayout,
  newton3,
  loadEstimator,
  traversal
  /* bestConfigurations(container,dataLayout,newton3,loadEstimator,traversal,nsMinRuntime,scenario) */
;

CREATE VIEW configWinners AS
SELECT
  container,
  dataLayout,
  newton3,
  loadEstimator,
  traversal,
  COUNT(*) AS Wins
FROM
  bestConfigurations
GROUP BY
  container,
  dataLayout,
  newton3,
  loadEstimator,
  traversal
ORDER BY
  COUNT(*) DESC
  /* configWinners(container,dataLayout,newton3,loadEstimator,traversal,Wins) */
;

CREATE VIEW measuredConfigs AS
SELECT DISTINCT
  container,
  traversal,
  dataLayout,
  newton3,
  loadEstimator
FROM
  Measurement;

CREATE VIEW loserConfigs AS
SELECT
  container,
  traversal,
  dataLayout,
  newton3,
  loadEstimator
FROM
  measuredConfigs
  NATURAL LEFT OUTER JOIN configWinners
WHERE
  Wins IS NULL
;

CREATE VIEW containerWinners AS
SELECT
  container,
  SUM(Wins)
FROM
  configWinners
GROUP BY
  container
ORDER BY
  SUM(Wins) DESC
  /* containerWinners(container,"SUM(Wins)") */
;

CREATE VIEW traversalWinners AS
SELECT
  traversal,
  SUM(Wins)
FROM
  configWinners
GROUP BY
  traversal
ORDER BY
  SUM(Wins) DESC
  /* traversalWinners(traversal, "SUM(Wins)") */
;

CREATE VIEW MeasurementWithRankInScenario AS
SELECT
  *,
  RANK() OVER(
    PARTITION BY scenario
    ORDER BY
      nsRuntime
  ) AS rank
FROM
  Measurement
ORDER BY
  rank
  /* MeasurementWithRankInScenario(scenario,container,cellSizeFactor,traversal,loadEstimator,dataLayout,newton3,iteration,nsRuntime,rank) */
;

CREATE VIEW MeasurementRankFactor AS
SELECT
  *,
  RANK() OVER sameScenarioRuntimeOrdered AS rank,
  CAST(nsRuntime AS DOUBLE) / MIN(nsRuntime) OVER sameScenarioRuntimeOrdered AS factorWorseThanBest
FROM
  Measurement WINDOW sameScenarioRuntimeOrdered AS (
    PARTITION BY scenario
    ORDER BY
      nsRuntime
  )
ORDER BY
  rank
  /* MeasurementRankFactor(scenario,container,cellSizeFactor,traversal,loadEstimator,dataLayout,newton3,iteration,nsRuntime,rank,factorWorseThanBest) */
;

CREATE VIEW configRanks AS
SELECT
  container,
  dataLayout,
  newton3,
  traversal,
  loadEstimator,
  MAX(rank) - MIN(rank) AS diffMinMaxRank,
  MAX(rank) AS maxRank,
  MIN(rank) AS minRank,
  AVG(rank) AS avgRank,
  MAX(factorWorseThanBest) / MIN(factorWorseThanBest) AS factorWorseToBestFactor,
  MAX(factorWorseThanBest) AS maxFactorWorseThanBest,
  MIN(factorWorseThanBest) AS minFactorWorseThanBest,
  AVG(factorWorseThanBest) AS avgFactorWorseThanBest
FROM
  MeasurementRankFactor
GROUP BY
  container,
  dataLayout,
  newton3,
  traversal,
  loadEstimator
  /* configRanks(container,dataLayout,newton3,traversal,loadEstimator,diffMinMaxRank,maxRank,minRank,avgRank) */
;

CREATE VIEW ScenarioTuningTime AS
SELECT
  scenario,
  SUM(nsRuntime) AS nsTuningTime
FROM
  Measurement
GROUP BY
  scenario
  /* ScenarioTuningTime(scenario,nsTuningTime) */
;

CREATE VIEW ScenarioWithDerived AS /* Replaces ScenarioTuningTime */
SELECT
  scenario,
  SUM(nsRuntime) AS nsTuningTime,
  MAX(rank) AS numConfigsMeasured,
  MIN(nsRuntime) AS nsMinIterationTime,
  SUM(nsRuntime) / CAST((MIN(nsRuntime) * MAX(rank)) AS DOUBLE) AS factorToOptimalRuntime
FROM
  MeasurementRankFactor
GROUP BY
  scenario
ORDER BY
  factorToOptimalRuntime
  /* ScenarioTuningTime(scenario,nsTuningTime,numConfigsMeasured,factorToOptimalRuntime) */
;

CREATE VIEW MeasurementRankFactorTuningTime AS
SELECT
  *,
  RANK() OVER sameScenarioRuntimeOrdered AS rank,
  CAST(nsRuntime AS DOUBLE) / MIN(nsRuntime) OVER sameScenarioRuntimeOrdered AS factorWorseThanBest,
  (
    SUM(CAST(nsRuntime AS DOUBLE)) OVER sameScenarioRuntimeOrdered
  ) / nsTuningTime AS percentTuningTimeUntil
FROM
  Measurement m
  NATURAL JOIN ScenarioTuningTime WINDOW sameScenarioRuntimeOrdered AS (
    PARTITION BY m.scenario
    ORDER BY
      nsRuntime
  )
ORDER BY
  rank
  /* MeasurementRankFactorTuningTime(scenario,container,cellSizeFactor,traversal,loadEstimator,dataLayout,newton3,iteration,nsRuntime,nsTuningTime,rank,factorWorseThanBest,percentTuningTimeUntil) */
;

CREATE VIEW PercentTuningTimePerRank AS
SELECT
  rank,
  MIN(percentTuningTimeUntil) AS minPercent,
  MAX(percentTuningTimeUntil) AS maxPercent,
  AVG(percentTuningTimeUntil) AS avgPercent
FROM
  MeasurementRankFactorTuningTime
WHERE
  ( SELECT
      numConfigsMeasured
    FROM
      ScenarioWithDerived s
    WHERE
      MeasurementRankFactorTuningTime.scenario = s.scenario
  ) = (SELECT * FROM numDifferentConfigs)
GROUP BY
  rank
ORDER BY
  rank
  /* PercentTuningTimePerRank(rank,minPercent,maxPercent,avgPercent) */
;

CREATE VIEW lastConfigsOfUnfinishedScenarios AS
SELECT
  *
FROM
  ScenarioWithDerived
  NATURAL JOIN MeasurementRankFactor
WHERE
  numConfigsMeasured < (SELECT * FROM numDifferentConfigs)
  AND numConfigsMeasured = rank
;

CREATE VIEW runtimeFactorBinning AS
SELECT
  CAST(maxFactorBin AS VARCHAR(10)) AS maxFactorBin,
  numConfigsUntil - LAG(numConfigsUntil, 1, 0) OVER (ORDER BY maxFactorBin ASC) AS configsInBucket
FROM (
	SELECT
	  maxFactorBin,
	  COUNT(*) AS numConfigsUntil
	FROM
		MeasurementRankFactor,
		(
			SELECT
			  column1 AS maxFactorBin
			FROM
			  (VALUES (1.0), (1.2), (1.5), (2), (3.5), (5), (7.5), (10), (15), (25), (50), (100), (250), (500), (1000), (2000), (1000000))
		)
	WHERE
		factorWorseThanBest <= maxFactorBin
	GROUP BY
		maxFactorBin
	)
;

CREATE VIEW uselessConfigsAndTheirSupercedingConfigs AS
SELECT
  worseConfig.container AS worseContainer,
  worseConfig.traversal AS worseTraversal,
  worseConfig.dataLayout AS worseDataLayout,
  worseConfig.newton3 AS worseNewton3,
  worseConfig.loadEstimator AS worseLoadEstimator,
  betterConfig.container AS betterContainer,
  betterConfig.traversal AS betterTraversal,
  betterConfig.dataLayout AS betterDataLayout,
  betterConfig.newton3 AS betterNewton3,
  betterConfig.loadEstimator AS betterLoadEstimator
FROM
  measuredConfigs betterConfig, measuredConfigs worseConfig
WHERE
  NOT EXISTS (
    SELECT
      *
    FROM
      MeasurementRankFactor worse,
      MeasurementRankFactor better
    WHERE
      worse.scenario = better.scenario

      AND worse.container = worseConfig.container
      AND worse.traversal = worseConfig.traversal
      AND worse.loadEstimator = worseConfig.loadEstimator
      AND worse.dataLayout = worseConfig.dataLayout
      AND worse.newton3 = worseConfig.newton3

      AND better.container = betterConfig.container
      AND better.traversal = betterConfig.traversal
      AND better.loadEstimator = betterConfig.loadEstimator
      AND better.dataLayout = betterConfig.dataLayout
      AND better.newton3 = betterConfig.newton3

      AND worse.rank <= better.rank
  )
;

)SQL";
};