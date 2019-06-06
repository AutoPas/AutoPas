# Contributing to AutoPas

**Thanks for contributing to AutoPas!** 

Please keep the following notes in mind while working.

## C++
####General Notes
* used cpp standard: C++14. If there is a peace of code, which could be done better using C++17, please add a comment `@todo C++17` including the C++17 version of the code.
* If there is no special reason, use smart pointers instead of raw pointers.
* Use AutoPas wrapper functions for OpenMP (src/autopas/utils/WrapOpenMP.h) instead of OpenMP functions to allow building without enabled OpenMP.
* Use `#pragma once` instead of header guards.
* Use `constexpr` instead of `#define`.

#### Code Formatting
Source code needs to be formatted in google code style using clang format version 6.0. Wrong formatted code is automatically rejected.

#### Comment Style
* All public methods and attributes need to be documented with doxygen comments.
* The first comment inside a comment block (`/** <comment> */`) is automatically treated as brief comment and needs to end with a dot. The `brief` keyword is omitted (please delete occurrences).
* All comments containing todos should be prefixed with `@todo`, thus they are visible in the global todo list.
* date format: dd.mm.yyyy (please replace occurrences of other formats)

## Commits
* Try to remember to format your commits properly according to the above mentioned rules. Its frustrating for everyone to have failing CI builds due to incorrect formatting.
* Please avoid using commits to save your unfinished work before switching branches, this pollutes the commit history. Please use `git stash` instead.

## Issues
Feel free to open issues for bugs, feature requests or optimization ideas.