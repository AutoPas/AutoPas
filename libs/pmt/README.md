# PMT: Power Measurement Toolkit

[![pipeline status](https://git.astron.nl/RD/pmt/badges/master/pipeline.svg)](https://git.astron.nl/RD/pmt/-/commits/master)
[![Latest Release](https://git.astron.nl/RD/pmt/-/badges/release.svg)](https://git.astron.nl/RD/pmt/-/releases)

PMT is a high-level software library capable of collecting power consumption
measurements on various hardware. The library provides a standard interface to
easily measure the energy use of devices such as CPUs and GPUs in critical
application sections.

# Installation

First clone the repository:

```
git clone --recursive https://git.astron.nl/RD/pmt.git
```

The `--recursive` flag makes sure that the git submodules are also cloned. If
you cloned the repository without this flag, you can initialize the submodules
as follows:

```
git submodule update --init
```

To build the software, run the following commands:

1. Set-up cmake in a build directory and cd into the directory, e.g.
   `~/pmt/build`
   - `cmake <source dir path> -DCMAKE_INSTALL_PREFIX=<install dir path>`
1. Optionally further configure cmake through interactive build settings or with
   command line variables
   - use `ccmake` and/or add `-DPMT_BUILD_<SENSOR>=<0 or 1>` to the `cmake`
     commandline to select which PMT are built.
1. make and install
   - `make install`

# Usage

Include the header file into your program, e.g.:

```
#include “pmt.h”
```

Depending on which PMT implementations you have selected during the build, you
can now initialize any PMT instance:

```
std::unique ptr<pmt::PMT> sensor(pmt::nvml::NVML::Create());
```

or use:

```
std::unique ptr<pmt::PMT> sensor(pmt::Create("nvml"));
```

Next, you can start measuring power using the common api as specified in
`pmt.h`:

```
pmt::State start, end;
start = sensor−>read();
...
end = sensor−>read () ;
std::cout<<sensor−>joules(start, end) <<” [J]“<<std::endl;
std::cout<<sensor−>watts(start, end) <<” [W]“<<std::endl;
std::cout<<sensor−>seconds(start, end)<<” [S]“<<std::endl;
```

## CMake

Integrating PMT into your CMake project is easy, as `.cmake` config files will
be automatically generated and put into your install directory during
installation. All you need to do is add `find_package(pmt)` to your
`CMakeLists.txt` and use `target_link_libraries(<target> pmt)` to setup the
include directories and link libraries. Add
`-Dpmt_ROOT=<pmt-installation-prefix>` to your CMake commandline to let CMake
find the PMT config files.

# PMT executable

The `PMT` executable might be used directly to read the values of a power meter
at a regular interval. The `PMT` executable is available in the install
directory `/bin`.

# Acknowledgement

If you decide to use PMT in your research, please cite the following reference:

```
@INPROCEEDINGS{10027520,
  author={Corda, Stefano and Veenboer, Bram and Tolley, Emma},
  booktitle={2022 IEEE/ACM International Workshop on HPC User Support Tools (HUST)},
  title={PMT: Power Measurement Toolkit},
  year={2022},
  volume={},
  number={},
  pages={44-47},
  doi={10.1109/HUST56722.2022.00011}
}
```
