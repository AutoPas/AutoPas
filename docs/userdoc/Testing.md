# Testing

AutoPas uses [googletest](https://github.com/google/googletest) as testing
framework and exposes tests to ctest, the CMake test driver.

## Running Tests

There are multiple possibilities. In order of recommendation:

1. Using `ctest`:
   ```bash
   ctest # add --verbose for more details on the tests
   ```
   To only run specific tests use arguments like `-R` (run tests matching regex) and `-E` (exclude tests matching regex)
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
   To only run specific tests use arguments
   ```bash
   tests/testAutopas/runTests --gtest_filter=ArrayMathTest.testAdd*
   ```

## Debugging Tests
Many IDEs (e.g., CLion) have integrated support for googletest, and you can debug the tests directly within the IDE.

If you prefer `gdb`:
1. Find out the command to start your desired test with `-N` aka. `--show-only`:
   ```bash
   ctest -R 'Array.*testAdd' -N
   ```
2. Start the test with `gdb`
   ```bash
   gdb --args ${TestCommand}
   ```