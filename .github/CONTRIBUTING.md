# Contributing to AutoPas

**Thanks for contributing to AutoPas!** 

Please keep in mind the following notes while working.

## C++
### General Notes
* Cpp standard: C++20. If there is a piece of code, which could be done better using a newer standard, please add a comment like `@todo C++23` including the alternative version of the code.
* Pointers: Always use smart pointers when you are managing memory. Don't use `new` or `delete`.
* OpenMP: Use AutoPas wrapper functions and macros for OpenMP from [`WrapOpenMP.h`](/src/autopas/utils/WrapOpenMP.h) instead of native OpenMP to allow clean building without OpenMP.
* `#pragma once` instead of header guards.
* `#include` of files from within AutoPas shall be given with the full path (starting with `autopas/`) and using `""`. 
* `constexpr` instead of `#define`. Use it wherever possible.
* `const` wherever possible. 
* `nullptr` instead of `NULL`.
* `using` instead of `typedef`.
* Avoid `assert()` but use `autopas::utils::ExceptionHandler::exception("Meaningful error message")` instead.

### Code Style
* Private attributes are prefixed with `_`.
* Postfix template class/typename parameters with `_T`. (This has not been applied retroactively and so there are plenty of cases where this is not yet the case, but should be applied for new contributions.
* Every (abstract) class gets its own file, named exactly like the class.
* Class names start with a capital letter.
* Use `camelCase` over `snake_case`.
* Google code style is enforced by the CI server.
* To enable code formatting targets set the `cmake` variable `AUTOPAS_FORMATTING_TARGETS` to `ON`.
* [Clang format](https://releases.llvm.org/14.0.0/tools/clang/docs/ClangFormat.html) version 14 is enforced (other versions might format slightly differently - you must use version 14 only).
* Run `make clangformat` before submitting a PR for review.
* [cmake format](https://github.com/cheshirekow/cmake_format) is enforced.
* Run `make cmakeformat` before submitting a PR for review.

### Comment Style
* Please write full sentences starting with a capital letter and ending with a period.
* [Doxygen](https://www.doxygen.nl/) (> v1.8.11) is used in this project to create the documentation.
* Documentation style is [Javadoc style](https://en.wikipedia.org/wiki/Javadoc).
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

## CI
A continuous integration setup (CI) is automatically run for each open pull request and for the master.
The executed tests are defined within the [GitHub workflow file](/.github/workflows/TestSuites.yaml).
These tests include:
* Formatting and documentation checks 
* Building of all targets and execution of the provided ctest tests.
* Sanitizer runs (Address+Leak sanitizer, Thread Sanitizer)

If you encounter problems within these tests check whether you can reproduce them locally. Have a look at the workflow file for how the tools are used with AutoPas. If you do not have the respective compiler installed you can use the [AutoPas docker images](https://hub.docker.com/u/autopas).
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
* Think about inheriting from a similar traversal. At least derive your new traversal from [`TraversalInterface`](/src/autopas/containers/TraversalInterface.h).
* Go to [`TraversalOption`](/src/autopas/options/TraversalOption.h).
  * Add a new enum in `TraversalOption::Value`.
  * Add a new string representation in the `map` of `TraversalOption::getOptionNames()`.
  * If the new traversal is a triwise traversal, add it to `TraversalOption::getAllTriwiseTraversals()`.
* If applicability of the traversal is restricted, add your new enum to any of the functions that return sets of restricted traversals in [`CompatibleTraversals`](/src/autopas/containers/CompatibleTraversals.h).
* Add a case for the new traversal in [`TraversalSelector::generateTraversal()`](/src/autopas/tuning/selectors/TraversalSelector.h).
* Check that the new option is working in the md-flexible example.
* Adapt unit tests (e.g. expected number of iterations in [`AutoTunerTest::testAllConfigurations()`](/tests/testAutopas/tests/tuning/AutoTunerTest.cpp)).
* Add new unit tests for your traversal.
* Regenerate the [`RuleLanguage.g4`](/src/autopas/tuning/tuningStrategy/ruleBasedTuning/RuleLanguage.g4) via the [`generateRuleLanguage.sh`](/src/autopas/tuning/tuningStrategy/ruleBasedTuning/generateRuleLanguage.sh) script, both located in [`ruleBasedTuning`](/src/autopas/tuning/tuningStrategy/ruleBasedTuning).

### Adding a new Container
* Create a new container class under [`src/autopas/containers/`](/src/autopas/containers/).
* Derive your new container from [`ParticleContainerInterface`](/src/autopas/containers/ParticleContainerInterface.h) or a more similar one.
* Go to [`ContainerOption`](/src/autopas/options/ContainerOption.h).
  * Add a new enum in `ContainerOption::Value`.
  * Add a new string representation in the `map` of `ContainerOption::getOptionNames()`.
* Create a new set of compatible traversals in [`CompatibleTraversals`](/src/autopas/containers/CompatibleTraversals.h).
* Create a new `case` statement in [`StaticContainerSelector`](/src/autopas/utils/StaticContainerSelector.h).
* Add a case for the new container in [`ContainerSelector::generateContainer()`](/src/autopas/tuning/selectors/ContainerSelector.h).
* Check that the new option is working in the md-flexible example.
* Adapt unit tests (e.g. expected number of iterations in [`AutoTunerTest::testAllConfigurations()`](/tests/testAutopas/tests/tuning/AutoTunerTest.cpp) and [`StringUtilsTest::parseContainerOptionsTest`](/tests/testAutopas/tests/utils/StringUtilsTest.cpp)).
* Add new unit tests for your container.
* Regenerate the [`RuleLanguage.g4`](/src/autopas/tuning/tuningStrategy/ruleBasedTuning/RuleLanguage.g4) via the [`generateRuleLanguage.sh`](/src/autopas/tuning/tuningStrategy/ruleBasedTuning/generateRuleLanguage.sh) script, both located in [`ruleBasedTuning`](/src/autopas/tuning/tuningStrategy/ruleBasedTuning).

### Adding a new Tuning Strategy
* Create a new tuning strategy class under [`src/autopas/tuning/tuningStrategy`](/src/autopas/tuning/tuningStrategy).
* Derive your new strategy from [`TuningStrategyInterface`](/src/autopas/tuning/tuningStrategy/TuningStrategyInterface.h) or a more similar one.
* Go to [`TuningStrategyOption`](/src/autopas/options/TuningStrategyOption.h).
  * Add a new enum in `TuningStrategyOption::Value`.
  * Add a new string representation in the `map` of `TuningStrategyOption::getOptionNames()`.
* In [`TuningStrategyFactory::generateTuningStrategy()`](/src/autopas/tuning/tuningStrategy/TuningStrategyFactory.cpp):
  * Add a `case` for the new strategy.
  * If the new strategy handles communication between processes itself, make sure to not wrap it in the lower `mpiStrategyOption`-switch.
* Check that the new option is working in the md-flexible example.
* Add new unit tests for your strategy.

### Adding a new Option
* If applicable add a new setter to [`AutoPas`](/src/autopas/AutoPas.h) (this is required for tunable options).
* Check that the new option is added to the md-flexible example. Parser and main.
* Global options, which are represented by an enum, should be defined in an additional file in [`src/autopas/options`](/src/autopas/options).
* Inherit from [`Option`](/src/autopas/options/Option.h). This will also generate functions for conversion from and to strings.
* Add new unit tests for your option, mainly in [`OptionTest`](/tests/testAutopas/tests/options/OptionTest.cpp).
* Also add the new option to md-flexible! 
  * Add the option in [`MDFlexConfig`](/examples/md-flexible/src/configuration/MDFlexConfig.h).
  * Add it to [`MDFlexConfig::to_string`](/examples/md-flexible/src/configuration/MDFlexConfig.h).
  * Parse it in [`CLIParser`](/examples/md-flexible/src/configuration/CLIParser.cpp).
  * In [`CLIParser::parseInput()`](/examples/md-flexible/src/configuration/CLIParser.cpp) add it to `relevantOptions` and switch.
  * Parse it in [`YamlParser`](/examples/md-flexible/src/configuration/YamlParser.cpp).
  * If applicable, pass the option value to AutoPas in [`Simulation::Simulation()`](/examples/md-flexible/src/Simulation.cpp).
  * Make sure that the description is parsable by [`CLIParser::createZSHCompletionFile()`](/examples/md-flexible/src/configuration/CLIParser.cpp).

### Making an Option tunable
* If not already done, add a new setter to [`src/autopas/AutoPasDecl.h`](/src/autopas/AutoPasDecl.h).
* Add your option to [`Configuration`](/src/autopas/tuning/Configuration.h) and adjust constructors, comparison operators and ConfigHash function accordingly.
* Adjust the individual tuning strategies accordingly; the exact implementation will depend on the purpose of your option, but some general advice is:
  * Depending on your new option, it might make sense for some tuning strategies to merge it with another option to avoid sparse dimensions.
  * Make sure it is added to the search spaces that is passed to the [`AutoTuner`](/src/autopas/tuning/AutoTuner.cpp) in `AutoPas::init()`
  * For Bayesian based tuning strategies your option will also have to be integrated into [`FeatureVector`](/src/autopas/tuning/utils/FeatureVector.h) and [`FeatureVectorEncoder`](/src/autopas/tuning/utils/FeatureVectorEncoder.h).
  * Extend [`FeatureVectorEncoder`](/src/autopas/tuning/utils/FeatureVectorEncoder.h) by modifying `setAllowedOptions()`, `convertToTunable()` and `convertFromTunable()`. If the new option wasn't merged with another one you may have to add a new index to `DiscreteIndices` or `ContinuousIndices`
  * Make sure to declare your option by calling `configureTuningParameter()` in [`ActiveHarmony::resetHarmony()`](/src/autopas/tuning/tuningStrategy/ActiveHarmony.cpp).
* In [`AutoPasConfigurationCommunicator`](/src/autopas/utils/AutoPasConfigurationCommunicator.h):
  * Change the size and (de-)serialization of SerializedConfiguration
  * Add the new option to all appropriate functions and adjust their functioning respectively.
* In [`ConfigurationAndRankIteratorHandler`](/src/autopas/utils/ConfigurationAndRankIteratorHandler.h):
  * Add the new option wherever appropriate.
  * If the new options depends on others, implement it similarly to traversals, containers, and load estimators.
* Adjust any tests that are affected by these changes. The following tests will definitely require changes:
  * [`AutoPasInterfaceTest`](/tests/testAutopas/tests/autopasInterface/AutoPasInterfaceTest.cpp)
  * [`AutoTunerTest`](/tests/testAutopas/tests/tuning/AutoTunerTest.cpp)
  * [`FeatureVectorTest`](/tests/testAutopas/tests/tuning/utils/FeatureVectorTest.cpp)
  * Tests for the individual tuning strategies. See files in [`tests/testAutopas/tests/tuning/tuningStrategy/`](/tests/testAutopas/tests/tuning/tuningStrategy/).

### Rule Based Tuning and Antlr
Whenever any option is added that is part of the tuning procedure, it has to be added to the Antlr parsing logic and grammar.
To update the grammar, simply run [`generateRuleLanguage.sh`](/src/autopas/tuning/tuningStrategy/ruleBasedTuning/generateRuleLanguage.sh) and overwrite the old [`RuleLangugage.g4`](/src/autopas/tuning/tuningStrategy/ruleBasedTuning/RuleLangugage.g4) file.
To update the parser, the easiest way to do this is via the CLion plugin for Antlr4. Just right-click the g4 file -> 'Configure ANTLR' (Output directory, input grammar file, namespace `AutopasGeneratedRuleSyntax`, language) and then right click the g4 again -> 'Generate ANTLR Recognizer'. Don't forget to apply clang-format afterward.

**WARNING** The CLion plugin and Antlr version must match!. For example for Antlr version 4.9.1 the plugin version 1.16 is needed, otherwise incompatible parser code is generated. See [Antlr plugin's GitHub page](https://github.com/antlr/intellij-plugin-v4/releases) for what is compatible and get the plugin's binary from the [jetbrains plugin webpage](https://plugins.jetbrains.com/plugin/7358-antlr-v4/versions).
