# Building AutoPas

## Requirements
* CMake 3.14 or newer
* a CMake generator target (`make` is tested)
* a C++20 compiler (gcc13, clang16, and ~~icpc 2019~~ are tested)
* OpenMP (comes with GCC, for Clang you need `libomp`)

Optional:
* For `tuningLogToSQL`: `libsqlite3`
* For rule based tuning and fuzzy tuning: `pkg-config`. By default, rule based tuning and fuzzy tuning are disabled but can be enabled via the CMake 
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

By default, AutoPas tunes for the best configuration according to runtime. For all Linux based systems, it is also possible to tune for the algorithm that consumes the least energy.
This is implemented via [PMT](https://git.astron.nl/RD/pmt) with customization carried out to remove performance overhead in energy measurements.
The changes with respect to parent repository can be found in `AutoPas/libs/patches/patch-file-pmt-for-autopas.patch`, and detail on the reason behind the changes can be found in the [bachelor's thesis](https://mediatum.ub.tum.de/doc/1760019/1760019.pdf) by Maximilian Praus.
In the current version, energy consumption can be measured via Intel's RAPL or with LIKWID.
RAPL is set as the default sensor and is forced to compile always.

To use energy tuning, energy measurements must be enabled using the CMake option:
```bash
cmake -DAUTOPAS_ENABLE_ENERGY_MEASUREMENTS=ON .. 
```

The RAPL sensor is forcibly enabled by AutoPas by forcing the cmake option `PMT_BUILD_RAPL` to `ON`. However, other available sensors like `LIKWID` can also be used by separately compiling them using the CMake option:
```bash
cmake -DPMT_BUILD_LIKWID=ON .. 
```
Note that user must have `LIKWID` installed before enabling this option.
Additionally, an energy sensor used for measurement needs to be selected during runtime by specifying the option `energy-sensor` in the input file. 
The default is set to `rapl` and will be used when no option is specified. 
However, when multiple sensors are compiled, user can choose any of the available sensors by specifying the `energy-sensor` option.

#### RAPL Permissions
Energy measurements using Intel RAPL require read access to the system powercap files located at:
```bash
/sys/class/powercap/intel-rapl*
```
By default, these files are readable only by the root user and simulation need to be run with `sudo`. This is not ideal.
For example, running MPI applications with `sudo` is strongly discouraged, the recommended workaround is to adjust the permissions of the RAPL files once prior to running the simulation:
```bash
sudo chmod -R a+r /sys/class/powercap/intel-rapl*
```
This command grants read-only access to all users, allowing the simulation to access RAPL energy counters without requiring sudo.
However, these permission changes apply only until the next system reboot.
After a reboot, the command must be executed again unless a persistent [udev](https://linuxconfig.org/tutorial-on-how-to-write-basic-udev-rules-in-linux) rule is configured.

**Note**: On most Linux clusters, such persistent configuration is typically not required, as system administrators already provide user-level access to RAPL counters.

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
