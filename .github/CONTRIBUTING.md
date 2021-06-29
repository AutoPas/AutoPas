# Contributing to AutoPas

**Thanks for contributing to AutoPas!** 

Please keep in mind the following notes while working.

## C++
### General Notes
* Cpp standard: C++17. If there is a piece of code, which could be done better using a newer standard, please add a comment like `@todo C++20` including the alternative version of the code.
* Pointers: Always use smart pointers when you are managing memory. Don't use `new` or `delete`.
* OpenMP: Use AutoPas wrapper functions for OpenMP (`src/autopas/utils/WrapOpenMP.h`) instead of OpenMP functions to allow building without enabled OpenMP.
* `#pragma once` instead of header guards.
* `#include` of files from within AutoPas shall be given with the full path (starting with `autopas/`) and using `""`. 
* `constexpr` instead of `#define`. Use it wherever possible.
* `const` wherever possible. 
* `nullptr` instead of `NULL`.
* `using` instead of `typedef`.
* Avoid `assert()` but use `autopas::utils::ExceptionHandler::exception("Descriptive error message")` instead.

### Code Style
* Private attributes are prefixed with `_`.
* Every (abstract) class gets its own file, named exactly like the class.
* Class names start with a capital letter.
* Use camelCase over snake_case.
* Google code style is enforced by the CI server.
* To enable code formatting targets set the `cmake` variable `AUTOPAS_FORMATTING_TARGETS` to `ON`.
* Clang format version 9 is enforced (Other versions might format slightly differently).
* Use `make clangformat` before submitting a PR.
* [cmake format](https://github.com/cheshirekow/cmake_format) is enforced.
* Use `make cmakeformat` before submitting a PR.

### Comment Style
* Please write full sentences starting with a capital letter and ending with a period.
* Doxygen (v1.8.11) is used in this project to create the documentation.
* Documentation style is Javadoc style.
* All public methods and attributes need to be documented.
* The first comment inside a comment block (`/** <comment> */`) is automatically treated as a brief comment and needs to end with a period. The `brief` keyword is omitted (please delete occurrences).
* ToDos: All comments containing todos should be prefixed with `@todo`, thus they are visible in the global todo list.
* Date format: dd.mm.yyyy (please replace occurrences of other formats as you come across them)

## GitHub
### Pull Requests
* If you want to contribute a new feature, resolve an issue, etc. please create a new branch or alternatively fork from `master` where you implement your solution.
* Create a pull request against `master` and fill out our Pull Request form.
* Please take your time to give a concise description of your solution and how you tested it.
* Link related PRs or issues.
* The CI server will check your PR for correct formatting and complete documentation. It also executes all unit tests and several examples with all supported compilers with and without sanitizers.
* To merge your PR back to master, at least one review of a project admin is required.

### Issues
* Feel free to open issues for bugs, feature requests or optimization ideas.
* Tag the issues you write appropriately.

### Commits
* Use meaningful commit messages.
* Please avoid using commits to save your unfinished work before switching branches, this pollutes the commit history. Please use `git stash` instead.

## Docker
You want to compile AutoPas with any sanitizers but do not have the appropriate compiler? Don't fret! There are docker containers for that. The containers are built from the [AutoPas-Dockerfiles repository](https://github.com/AutoPas/AutoPas-Dockerfiles) and prebuilds are hosted at [dockerhub](https://hub.docker.com/search?q=autopas%2F&type=image). To use a compiler from a container either mount your AutoPas folder and start bash in the container:
```bash
docker run -v ${PathToAutoPasRoot}/:/autopas -it autopas/autopas-build-intel bash
```
or directly start the compilation process:
```bash
docker run -v ${PathToAutoPasRoot}/:/autopas -it autopas/autopas-build-intel \
  bash -c "cd /autopas/build \
  && cmake -G Ninja .. \
  && ninja"
```
The tests executed through Jenkins are using these docker images.

## Jenkins/CI
A continuous integration setup (CI) is automatically run for each open pull request and for the master.
The executed tests are defined within the Jenkinsfile in the root AutoPas directory.
These tests include:
* Formatting and documentation checks 
* Builing of all targets and execution of the provided ctest tests.
* Sanitizer runs (Address+Leak sanitizer, Thread Sanitizer)

If you encounter problems within these tests check whether you can reproduce them locally. Have a look at the [Jenkinsfile](https://github.com/AutoPas/AutoPas/blob/master/Jenkinsfile) for how the tools are used with AutoPas. If you do not have the respective compiler installed you can use the [AutoPas docker images](https://hub.docker.com/u/autopas).
To circumvent "unknown module" problems with the thread sanitizer, a library to override `dlclose()` can be found in the `libs` directory. Said library can be used to get better stack traces that are caused by dynamically loaded libraries (using `dlopen()`).
More details can be found [here](../libs/fake-dlclose/README.md).

## AutoPas

### Guidelines
* Try to require as few things from the `Particle` classes and from functors as possible.
* This includes that we do not make restrictions on the constructors of the Particle class, this is tested with `DifferentParticlesTest`. 

### Namespaces
* Code in folder `src` should belong to namespace `autopas`.
* Classes which shouldn't be used externally should belong to namespace `internal`.

### Logging
AutoPas has its own logger based on [spdlog](https://github.com/gabime/spdlog) which can be used after the initialization of an AutoPas object via:
```C++
AutoPasLog(warn, "Hello {}", name);
```
The global log level can be set at runtime with:
```C++
#include "autopas/utils/Logger.h"
autopas::Logger::get()->set_level(autopas::Logger::LogLevel::debug);
```
Possible log levels are:`trace`, `debug`, `info`, `warn`, `err`, `critical`, `off`,

### Adding a new Traversal
* Create a new traversal class under `src/autopas/containers/[Container]/traversals` for the container the traversal is intended.
* Think about inheriting from a similar traversal. At least derive your new traversal from `src/autopas/containers/cellPairTraversals/TraversalInterface.h`.
* Go to `src/autopas/options/TraversalOption.h`.
  * Add a new enum in `TraversalOption::Value`.
  * Add a new string representation in the `map` of `TraversalOption::getOptionNames()`.
* Add the enum to every compatible container in `src/autopas/containers/CompatibleTraversals.h`.
* Add a case for the new traversal in `src/autopas/selectors/TraversalSelector.h::generateTraversal()`.
* Check that the new option is working in the md-flexible example.
* Adapt unit tests (e.g. expected number of iterations in `tests/testAutopas/tests/selectors/AutoTunerTest.cpp::testAllConfigurations()` and `OptionTest::parseTraversalOptionsTest`).
* Add new unit tests for your traversal.

### Adding a new Container
* Create a new container class under `src/autopas/containers/`.
* Derive your new container from `src/autopas/containers/ParticleContainerInterface.h` or a more similar one.
* Go to `src/autopas/options/ContainerOption.h`.
  * Add a new enum in `ContainerOption::Value`.
  * Add a new string representation in the `map` of `ContainerOption::getOptionNames()`.
* Create a new set of compatible traversals in `src/autopas/containers/CompatibleTraversals.h`.
* Create a new `case` statement in `src/autopas/utils/StaticContainerSelector.h`.
* Add a case for the new container in `src/autopas/selectors/ContainerSelector.h::generateContainer()`.
* Check that the new option is working in the md-flexible example.
* Adapt unit tests (e.g. expected number of iterations in `tests/testAutopas/tests/selectors/AutoTunerTest.cpp::testAllConfigurations()` and `StringUtilsTest::parseContainerOptionsTest`).
* Add new unit tests for your container.

### Adding a new Tuning Strategy
* Create a new tuning strategy class under `src/autopas/selectors/tuningStrategy`.
* Derive your new strategy from `src/autopas/selectors/tuningStrategy/TuningStrategyInterface.h` or a more similar one.
* Go to `src/autopas/options/TuningStrategyOption.h`.
  * Add a new enum in `TuningStrategyOption::Value`.
  * Add a new string representation in the `map` of `TuningStrategyOption::getOptionNames()`.
* In `src/autopas/AutoPas.h::generateTuningStrategy()`:
  * Add a `case` for the new strategy.
  * If the new strategy handles communication between processes itself, make sure to not wrap it in the lower `mpiStrategyOption`-switch.
* Check that the new option is working in the md-flexible example.
* Add new unit tests for your strategy.

### Adding a new Option
* If applicable add a new setter to `src/autopas/AutoPas.h` (this is required for tunable options).
* Check that the new option is added to the md-flexible example. Parser and main.
* Global options, which are represented by an enum, should be defined in an additional file in `src/autopas/options`.
* Inherit from `src/autopas/options/Option.h`. This will also generate functions for conversion from and to strings.
* Add new unit tests for your option, mainly in `tests/testAutopas/tests/options/OptionTest.cpp`.
* Also add the new option to md-flexible! 
  * The option needs to be added in `examples/md-flexible/src/parsing/MDFlexConfig.h`.
  * Parsing for it in `examples/md-flexible/src/parsing/CLIParser.cpp` and `YamlParser.cpp`
  * Make sure that the description is parsable by `CLIParser::createZSHCompletionFile()`

### Making an Option tunable
* If not already done, add a new setter to `src/autopas/AutoPas.h`.
* Initiate the set of allowed options to all options in the constructor of AutoPas.
* Add your option to `src/autopas/selectors/Configuration.h` and adjust constructors, comparison operators and ConfigHash function accordingly.
* Add a parameter for your option to `TuningStrategyFactory::generateTuningStrategy` and pass it to the constructors for each tuning strategy.
* Adjust the individual tuning strategies accordingly; the exact implementation will depend on the purpose of your option, but some general advice is:
  * Depending on your new option, it might make sense for some tuning strategies to merge it with another option to avoid sparse dimensions.
  * `FullSearch` and `PredictiveTuning` inherit from `SetSearchSpaceBasedTuningStrategy`, adjust the constructor for this class and the method `populateSearchSpace()`.
  * For bayesian based tuning strategies your option will also have to be integrated into `FeatureVector` and `FeatureVectorEncoder`.
  * Extend `FeatureVectorEncoder` by modifying `setAllowedOptions()`, `convertToTunable()` and `convertFromTunable()`. If the new option wasn't merged with another one you may have to add a new index to `DiscreteIndices` or `ContinuousIndices`
  * Make sure to declare your option by calling `configureTuningParameter()` in `ActiveHarmony::resetHarmony()`.
* In `src/autopas/utils/AutoPasConfigurationCommunicator.h/.cpp`:
  * Change the size and (de-)serialization of SerializedConfiguration
  * Add the new option to all appropriate functions and adjust their functioning respectively.
* In `src/autopas/utils/ConfigurationAndRankIteratorHandler.h/.cpp`:
  * Add the new option wherever appropriate.
  * If the new options depends on others, implement it similarly to traversals, containers, and load estimators.
* Adjust any tests that are affected by these changes. The following tests will definitely require changes:
  * `tests/testAutopas/tests/autopasInterface/AutoPasInterfaceTest.{h,cpp}`
  * `tests/testAutopas/tests/selectors/AutoTunerTest.cpp`
  * `tests/testAutopas/tests/selectors/FeatureVectorTest.cpp`
  * Tests for the individual tuning strategies. See files in `tests/testAutopas/tests/selectors/tuningStrategy/`.
