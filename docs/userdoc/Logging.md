# Logging

AutoPas contains multiple loggers with different purposes that can help to shed light into the black box.
Under the hood, they use [spdlog](https://github.com/gabime/spdlog).
When deactivated via `CMake` these loggers do not add any run time overhead.

## AutoPasLog
This is the main, general purpose logger. It supports all spdlog-levels. These levels can be (de)activated at compile time via the `CMake` variable `AUTOPAS_MIN_LOG_LVL`. At run time, this logger's compiled levels can be set e.g. via:
`autopas::Logger::get()->set_level(autopas::Logger::LogLevel::debug);`

At debug level, this logger will print the full configuration of every call to `autopas::AutoPas::iteratePairwise()`.

## GaussianClusterLogger
Creates a graph representation of the Gaussian cluster model that was created during the simulation.
This logger is switched on/off via the `CMake` variable `AUTOPAS_LOG_GAUSSIANCLUSTER`.

## IterationLogger
Creates a csv file containing information about the configuration and timings of every single pairwise iteration.
This logger is switched on/off via the `CMake` variable `AUTOPAS_LOG_ITERATIONS`.

## OctreeLogger
Creates a vtk file to visualize the Octree particle container.
This logger is switched on/off via the `CMake` variable `AUTOPAS_LOG_OCTREE`.

## PredictionLogger
Creates a csv containing the predictions made by the PredictiveTuning strategy.
This logger is switched on/off via the `CMake` variable `AUTOPAS_LOG_PREDICTIONS`.

## TuningDataLogger
Creates a csv containing all data that is collected for tuning purposes. This is the raw data that is available to
the tuning algorithms.
This logger is switched on/off via the `CMake` variable `AUTOPAS_LOG_TUNINGDATA`.

## TuningResultLogger
Creates a csv containing the results of every tuning phase. Useful if only the high level end results are of interest.
This logger is switched on/off via the `CMake` variable `AUTOPAS_LOG_TUNINGRESULTS`.
