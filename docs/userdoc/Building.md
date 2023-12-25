# Building AutoPas

## Requirements
* CMake 3.14 or newer
* make (build-essentials) or ninja
* a C++17 compiler (gcc11, clang13, and ~~icpc 2019~~ are tested.)
* For rule based tuning: `pkg-config`, `uuid`
* For `tuningLogToSQL`: `libsqlite3`

## Build Instructions
build instructions for make:
```bash
mkdir build
cd build
cmake ..
make
```
if you want to use a specific compiler, specify it at the first CMake call, e.g.:
```bash
mkdir build
cd build
CC=clang CXX=clang++ cmake ..
make
```
if you would like to use ninja instead of make:
```bash
mkdir build
cd build
cmake -G Ninja ..
ninja
```
### Building AutoPas on a Cluster
HPC clusters often use module systems. CMake is sometimes not able to
correctly detect the compiler you wished to use. If a wrong compiler is
found please specify the compiler explicitly, e.g. for gcc:
```bash
mkdir build
cd build
CC=`which gcc` CXX=`which g++` cmake ..
make
```

### Dependency Management
AutoPas relies on a small number of dependencies. By default, AutoPas looks for
installed versions of those libraries, but it can also be forced to (selectively)
use bundled versions. To make use of this feature, call `cmake` with:
```bash
cmake -D spdlog_ForceBundled=ON    # replace spdlog by the lib you want to force
```
For a full list, have a look at the variables exposed via `ccmake`. 

## Profiling Compilation Time 

For profiling the compile-time, the `cmake` option `AUTOPAS_COMPILE_TIME_PROFILING` can be turned on. This enables gcc's -`ftime-report` and clang's `-ftime-trace`.
It is recommended to use clang, as its output is more detailed.
`-ftime-trace` generates a .json file for each compilation unit next to the generated object file (inside one of the CMakeFiles directories).
Chrome has a built-in tool for viewing these files in a flame graph. It can be accessed through the URL `chrome://tracing`.
