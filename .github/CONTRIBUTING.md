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
* Avoid `assert()` but use `autopas::utils::ExceptionHandler::exception("Descriptive error message")` instead.

### Code Style
* Private attributes are prefixed with `_`.
* Every (abstract) class gets its own file, named exactly like the class.
* Class names start with a capital letter.
* Use camelCase over snake_case.
* Google code style is enforced by the CI server.
* Clang format version 6.0 (Other versions might format slightly differently).
* Use `make clangformat` before submitting a PR.

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

## AutoPas
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
* Add a new enum entry to `src/autopas/options/TraversalOption.h`.
* Add the enum to every compatible container in `src/autopas/containers/CompatibleTraversals.h`.
* Add new parsing and toString cases to `src/autopas/utils/StringUtils.h`.
* Add a case for the new traversal in `src/autopas/selectors/TraversalSelector.h::generateTraversal()`.
* Check that the new option is added to the md-flexible example.
* Adapt unit tests (e.g. expected number of iterations in `tests/testAutopas/tests/selectors/AutoTunerTest.cpp::testAllConfigurations()` and `StringUtilsTest::parseTraversalOptionsTest`).
* Add new unit tests for your traversal.

### Adding a new Container
* Create a new container class under `src/autopas/containers/`.
* Derive your new container from `src/autopas/containers/ParticleContainer.h` or a more similar one.
* Add a new enum entry to `src/autopas/options/ContainerOption.h`.
* Create a new set of compatible traversals in `src/autopas/containers/CompatibleTraversals.h`.
* Create a new `case` statement in `src/autopas/utils/StaticSelectors.h`.
* Add new parsing and toString cases to `src/autopas/utils/StringUtils.h`.
* Add a case for the new container in `src/autopas/selectors/ContainerSelector.h::generateContainer()`.
* Check that the new option is added to the md-flexible example.
* Adapt unit tests (e.g. expected number of iterations in `tests/testAutopas/tests/selectors/AutoTunerTest.cpp::testAllConfigurations()` and `StringUtilsTest::parseContainerOptionsTest`).
* Add new unit tests for your container.

### Adding a new Tuning Strategy
* Create a new tuning strategy class under `src/autopas/selectors/tuningStrategy`.
* Derive your new strategy from `src/autopas/selectors/tuningStrategy/TuningStrategyInterface.h` or a more similar one.
* Add a new enum entry to `src/autopas/options/TuningStrategyOption.h` along with a brief description.
* Add new parsing and toString cases to `src/autopas/utils/StringUtils.h`.
* Add a case for the new strategy in `src/autopas/AutoPas.h::generateTuningStrategy()`.
* Check that the new option is added to the md-flexible example.
* Add new unit tests for your strategy.

### Adding a new Option
* If applicable add a new setter to `src/autopas/AutoPas.h`.
* Check that the new option is added to the md-flexible example. Parser and main.
* Global options, which are represented by an enum, should be defined in an additional file in `src/autopas/options`.
* Add new unit tests for your option.
