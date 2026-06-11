# Developer documentation

If you're looking for user documentation, go [here](README.md).

Before you can do development work, you'll need to check out a local copy of the
repository:

```shell
cd <where you keep your git repositories>
git clone https://git.astron.nl/RD/pmt.git
cd pmt
```

## Prerequisites

### Development environment

Any modern Linux system should suffice as a development environment.

### Build tools

Summary of what you need :

- `gcc` 9 or above
- `g++` 9 or above
- `make` 4 or above
- `cmake` 3.20 or above

Check that you have the correct `gcc`, `g++` and `make` versions using

```shell
gcc --version
g++ --version
make --version
cmake --version
```

On a Debian-based system you can install them with

```shell
sudo apt install build-essential
```

Next, you need CMake 3.20 or above. Check if you have the correct version
installed with `cmake --version`. If your CMake version is not adequate, you can
install CMake manually by downloading the latest **stable** version from the
[CMake downloads page](https://cmake.org/download/) and following the
[installation instructions](https://cmake.org/install/).

If you don't have enough privileges to install `cmake` globally - for instance
if you are in a cluster without privileges - you can use `--prefix=PREFIX` to
install the CMake to your home folder. Remember that your `PATH` variable must
contain the path that you install `cmake`, so for instance, you can add the
following to your `.bashrc`:

```shell
PREFIX=<PREFIX-USED-WITH-CMAKE>
export PATH=$PREFIX/bin:$PATH
```

Remember to update your environment either by logging out and in again, or
running `source $HOME/.bashrc`.

## Linting and Formatting

We use the following linters and formatters in this project:

- [clang-format](https://clang.llvm.org/docs/ClangFormat.html)
- [cmake-format](https://cmake-format.readthedocs.io/en/latest/installation.html)
- [include-what-you-use](https://github.com/include-what-you-use/include-what-you-use)

The formatter `clang-format` will format all source files, and `cmake-format`
will format all CMake-related files. `include-what-you-use` will check whether
you included all the required header files.

To run the formatters and linters, you first need to build the project and set
up \`pre-commit\`\`. See [this](#pre-commit-hook) section for more details.

### pre-commit hook

`pre-commit` is a tool that can automatically run linters, formatters, or any
other executables whenever you commit code with `git commit`.

If you think having such automated checks is helpful for development, you can
install the pre-commit CLI from PyPI using pip:

```shell
python3 -m pip install --user pre-commit
```

For other install options, look [here](https://pre-commit.com/#installation).

You can install formatting and linting tools as follows:

```shell
python3 -m pip install --user clang-format cmake-format
apt install iwyu
```

Enable the pre-commit hooks defined in `.pre-commit-config.yaml` with:

```shell
pre-commit install
```

Once enabled, future `git commit`s will trigger the pre-commit hooks. Depending
on which files are changed by a given commit, some checks will be skipped.

You can uninstall the pre-commit hooks by:

```shell
pre-commit uninstall
```

The `include-what-you-use` extension needs the `compile_commands.json` from the
build directory. Make sure that there is `build` folder in your source and tree
and that you first ran `cmake`. `build` can also be a symbolic link for
out-of-tree builds. Thus if you cloned pmt into `~/src/pmt` and use
`~/build/pmt` as a build directory, run:

```shell
ln -s ~/build/pmt ~/src/pmt/build
```

Running pre-commit hooks individually is also possible with:

```shell
pre-commit run <name of the task>
```

For example:

```shell
pre-commit run clang-format
pre-commit run cmake-format
pre-commit run include-what-you-use
```

To run all the tools at once:

```shell
pre-commit run -a
```
