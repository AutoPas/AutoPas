# Contributing to AutoPas

**Thanks for contributing to AutoPas!** 

Please keep in mind the following notes while working.

## C++
### General Notes
* Cpp standard: C++14. If there is a piece of code, which could be done better using a newer standard, please add a comment like `@todo C++17` including the alternative version of the code.
For now we are stuck with C++14 due to the CUDA dependency.
* Pointers: Unless there is no special reason, use smart pointers instead of raw pointers.
* OpenMP: Use AutoPas wrapper functions for OpenMP (src/autopas/utils/WrapOpenMP.h) instead of OpenMP functions to allow building without enabled OpenMP.
* `#pragma once` instead of header guards.
* `constexpr` instead of `#define`.

### Code Formatting
* Google code style is enforced by the CI server.
* Clang format version 6.0 (Other versions might format slightly differently).
* Use `make clangformat` before submitting a PR.

### Comment Style
* Doxygen is used in this project to create a documentation.
* Documentation style is Javadoc style.
* All public methods and attributes need to be documented.
* The first comment inside a comment block (`/** <comment> */`) is automatically treated as brief comment and needs to end with a dot. The `brief` keyword is omitted (please delete occurrences).
* ToDos: All comments containing TpDos should be prefixed with `@todo`, thus they are visible in the global todo list.
* Date format: dd.mm.yyyy (please replace occurrences of other formats as you come accross them)

## Commits
* Use meaningfull commit messages.
* Please avoid using commits to save your unfinished work before switching branches, this pollutes the commit history. Please use `git stash` instead.

## Issues
* Feel free to open issues for bugs, feature requests or optimization ideas.
* Tag the issues you write appropriately.
