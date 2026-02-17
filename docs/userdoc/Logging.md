# Logging

AutoPas contains multiple loggers with different purposes that can help shed light on the black box.
Under the hood, they use [spdlog](https://github.com/gabime/spdlog).
When deactivated via `CMake`, these loggers do not add any run time overhead.

## AutoPasLog
This is the main general-purpose logger.
It supports all spdlog-levels, which are also accessible via the type alias `autopas::Logger::LogLevel`.
These levels can be (de)activated at compile time via the `CMake` variable `AUTOPAS_MIN_LOG_LVL`.
At run time, this logger's compiled levels can be set, e.g., via:

```c++
autopas::Logger::get()->set_level(autopas::Logger::LogLevel::debug);
```

To conveniently use this logger, a macro is provided that takes the log level in all caps, a message string which will be formatted by [fmt](https://github.com/fmtlib/fmt) followed by any number of arguments:
```c++
AutoPasLog(INFO, "Message with two arguments: {} {}", arg0, arg1);
```

## GaussianClusterLogger
Creates a graph representation of the Gaussian cluster model that was created during the simulation.
This logger is switched on/off via the `CMake` variable `AUTOPAS_LOG_GAUSSIANCLUSTER`.

## IterationLogger
Creates a CSV file containing information about the configuration and timings of every single call to `computeInteractions()`.
This logger is switched on/off via the `CMake` variable `AUTOPAS_LOG_ITERATIONS`.

## OctreeLogger
Creates a VTK file to visualize the Octree particle container.
This logger is switched on/off via the `CMake` variable `AUTOPAS_LOG_OCTREE`.

## PredictionLogger
Creates a CSV containing the predictions made by the PredictiveTuning strategy.
This logger is switched on/off via the `CMake` variable `AUTOPAS_LOG_PREDICTIONS`.

## TuningDataLogger
Creates a CSV containing all data that is collected for tuning purposes.
This is the raw data that is available to
the tuning algorithms.
This logger is switched on/off via the `CMake` variable `AUTOPAS_LOG_TUNINGDATA`.

## TuningResultLogger
Creates a CSV containing the results of every tuning phase.
It is beneficial if only the high-level end results are of interest.
This logger is switched on/off via the `CMake` variable `AUTOPAS_LOG_TUNINGRESULTS`.

## FLOPLogger
Creates a CSV containing the "useful" FLOP count and hit rate per particle interaction.
By "useful" FLOPs, we mean avoiding anything which does not contribute to the particle interaction 
e.g. FLOPs spent on masked vector registers.
As this is subjective, please refer to the documentation of the Functor you are using for more details.
By hit rate, we mean the proportion of distance calculations that are within the cutoff.
If a functor has not implemented `getNumFLOPs` or `getHitRate`, "Not Implemented" will be outputted instead.
This logger is switched on/off via the `CMake` variable `AUTOPAS_LOG_FLOPS`.
