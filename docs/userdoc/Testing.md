# Testing

AutoPas uses [googletest](https://github.com/google/googletest) as its testing framework, which exposes tests to [`ctest`](https://cmake.org/cmake/help/latest/manual/ctest.1.html), the CMake test driver.

We do not distinguish between unit/integration/smoke/etc. in the tests directory `/tests/`. However, our tests, in general, fall into a few categories:
- Regular GTests

  Everything in [`/tests/testAutopas/`](https://github.com/AutoPas/AutoPas/blob/master/tests/testAutopas/)

- GTests with MPI

  Everything in [`/tests/MPIParallelAutoPasTests/`](https://github.com/AutoPas/AutoPas/blob/master/tests/MPIParallelAutoPasTests/)

- Test involving examples

  Declared in the examples' CMakeLists files. [E.g. md-flexible](https://github.com/AutoPas/AutoPas/blob/master/examples/md-flexible/CMakeLists.txt).

- GTests for examples
  This only exists for md-flexible in `examples/md-flexible/tests`. See md-flexible's [README](https://github.com/AutoPas/AutoPas/blob/master/examples/md-flexible/README.md) for more information.

## Running Tests

There are multiple possibilities to run tests from the `/tests/` directory. In order of recommendation:

1. Using `ctest`:
   ```bash
   ctest # add --verbose or --output-on-failure for more details on the tests and -j N to run N tests in parallel 
   ```
   To only run specific tests, use arguments like `-R` (run tests matching regex) and `-E` (exclude tests matching regex)
   ```bash
   ctest -R 'Array.*testAdd' -E 'Double'
   ```
2. Using the `make` target:
   ```bash
   make test
   ```
3. Directly launching the test executable:
   ```bash
   tests/testAutopas/runTests
   ```
   To only run specific tests, use arguments
   ```bash
   tests/testAutopas/runTests --gtest_filter=ArrayMathTest.testAdd*
   ```

### Running Example Tests

To enable the example simulations, run
```bash
ctest -C checkExamples
```
Example tests defined in `CMake` files are disabled by default.
If you only want to run only those specific to one example, in your build directory, decent into its directory and run the command from above there.

## Debugging Tests
Many IDEs (e.g., CLion) have integrated support for googletest, and you can debug the tests directly within the IDE.

If you prefer `gdb` from the command line:
1. Find the command to start your desired test with `-N` (aka. `--show-only`) plus `-V` (aka. `--verbose`):
   ```bash
   ctest -R 'Array.*testAdd' -N -V
   ```
2. Start the test with `gdb`
   ```bash
   gdb --args ${TestCommand}
   ```

### MPI tests

MPI parallel applications like md-flexible can be a bit tricky to debug, so make sure you debug without MPI whenever you can.
If you have to debug MPI code, the trick is to attach one `gdb` per rank and use as few ranks as needed.
Remember to also limit the number of OpenMP threads to not overload your system.

```bash
# Maximum number of threads for an eight-core processor
OMP_NUM_THREADS=2 mpirun -n 4 -xterm -1! gdb --args examples/md-flexible/md-flexible
```
Explanation:

| Command Snippet     | Description                                                                                                                                                                                                                        |
|:-------------------:|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `OMP_NUM_THREADS=2` | Shell environment variable to lock OpenMP to threads per rank (=process).                                                                                                                                                          |
| `mpirun -n 4`       | Four ranks will be started.                                                                                                                                                                                                        |
| `-xterm -1!`        | Open an `xTerm` console window for ranks. `-1` indicates that it should be done for all ranks, no matter the count, `!` ensures the windows do not immediately close when the program stops or crashes so that you can read the output. |
| `gdb --args`        | Launches the debugger. Everything that follows is treated as arguments to `gdb`.                                                                                                                                                   |

## Creating Coverage Reports
The following commands can be used to create a coverage report in the form of HTML output. Please note that this is only supported with GCC, and debugging mode must be enabled to create coverage reports. Furthermore, LCOV version >= 2.0 must be installed.

```bash
cd $BUILD_DIR
cmake -DCMAKE_BUILD_TYPE=Debug -DAUTOPAS_BUILD_TESTS=ON -DAUTOPAS_ENABLE_COVERAGE=ON ..
make coverage
```

After these commands have been executed, a folder named `$BUILD_DIR/coverage` can be found which contains the HTML output. This can be downloaded and viewed in a browser.

## Related Files and Folders