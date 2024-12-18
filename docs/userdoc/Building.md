# Building AutoPas

## Requirements
* CMake 3.14 or newer
* a CMake generator target (`make` is tested)
* a C++17 compiler (gcc11, clang13, and ~~icpc 2019~~ are tested)
* OpenMP (comes with GCC, for Clang you need `libomp`)

Optional:
* For `tuningLogToSQL`: `libsqlite3`
* For rule based tuning and fuzzy tuning: `automake` (for `uuid`), `pkg-config`. By default, rule based tuning and fuzzy tuning are disabled but can be enabled via the CMake 
(see [below](#rules-based-tuning-fuzzy-tuning))

There are a few more dependencies, however you don't need to install them because they come bundled with AutoPas.
See [libs/](/libs) for a complete list.

If you insist on using your locally installed version of any of them, set the CMake variable `<LIBNAME>_ForceBundled=OFF`

## Build Instructions
Create a build directory and run the default CMake-based build workflow:
```bash
mkdir build
cd build
cmake ..    # Better: use the interactive version ccmake
cmake --build .
```

### Faster Compilation
It is highly recommended to only build the target you need instead of all default targets. E.g.
The following CMake features are highly recommended to keep compile time and system load to the necessary minimum:
- Only build the target you need
- Use multiple threads for compilation.
  How many threads you should use depends on how many threads your CPU supports and how much RAM you have available.
- Clang tends to use significantly less RAM per thread, thus enabling more parallelism

Here is an example of a parallel compilation of the md-flexible example:
```bash
cmake --build . --target md-flexible --parallel 12
```

### Enabling Rules-Based Tuning and Fuzzy Tuning
<a id="rules-based-tuning-fuzzy-tuning"></a>


Two of the possible tuning strategies of AutoPas, rules-based tuning and fuzzy-tuning, require dependencies `antlr4cpp` and `uuid`. These
are bundled with AutoPas, but can take some time to compile and, in some rare cases, lead to compilation errors. 
As such, by default, rule based tuning is disabled, and `antlr4cpp` and `uuid` are not compiled.

Both tuning strategies can be enabled via the CMake option:
```bash
cmake -DAUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING=ON .. 
```

### Energy Measurements and Tuning

By default, AutoPas tunes for the best configuration according to runtime. For all Linux based systems, it is also possible to tune
for the algorithm that consumes the least energy. This is implemented via [Intel's RAPL](https://www.intel.com/content/www/us/en/developer/articles/technical/software-security-guidance/advisory-guidance/running-average-power-limit-energy-reporting.html) 
(Running Average Power Limit) interface.

To use energy tuning, energy measurements must be enabled using the CMake option:
```bash
cmake -DAUTOPAS_ENABLE_ENERGY_MEASUREMENTS=ON .. 
```

### Select a Non-Default Compiler
If you want to use a different compiler than your system default, change the `CC` and `CXX` environment variables during initial configuration AND building:
```bash
CC=clang CXX=clang++ cmake ..
CC=clang CXX=clang++ cmake --build .
```

This assumes `clang` and `clang++` is in your search path.
Otherwise, specify the full path e.g. `/usr/bin/clang` or wherever your clang is kept.

Explicitly selecting the compiler is very often especially important when compiling on a cluster!

## Profiling Compilation Time 

AutoPas offers easy access to compiler-based compile-time profiling for GCC and Clang via the CMake option `AUTOPAS_COMPILE_TIME_PROFILING`.
This generates one `.json` file per compilation unit next to the generated object file deep in the build folder.
Chrome has a built-in tool for viewing these files as flame graphs.
It can be accessed through the URL `chrome://tracing`.

## Related Files and Folders
- cmake/
- CMakeLists.txt
