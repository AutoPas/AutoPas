# Power Measurement Toolkit (PMT) Releases

This repository contains releases of the Power Measurement Toolkit (PMT), a
comprehensive toolkit for power measurement and monitoring. Below are details
for each release:

## Unreleased

### Changed:

- Set NVML device by UUID

## 1.3.1

### Changed:

- Bugfixes
- Use instantaneous power measurement in `NVML`
- Change output of polling mode for `PMT` to show individual measurements

### Added:

- `NVML` measures GPU and module seperately (for Grace Hopper)
- `State::NrMeasurements()` and `State::{timetamp,name,joules,watts}(int)`

## 1.3.0

### Added:

- `pmt::Create` interface that accepts the platform as a string
- `pmt::NVIDIA` that automatically dispatches either `NVML` or `Tegra`
- Support to include PMT in another CMake project using `FetchContent` (see
  `example`)
- Add Likwid support based on performance groups (thanks @markbuettner!)

### Changed:

- The CMake options are prefixed with `PMT_`, rather than having it as `_PMT`
  suffix
- Simpler CMake code for `ROCM`
- The Python interface uses the new `pmt::Create` interface

## 1.2.0

### Changelog:

- Major refactoring + cleanup throughout (both in the library and in CMake)
- Add `PowerSensor3`
- Add `Cray`
- Remove deprecated `AMDGPU` (use `ROCM` instead)
- Reimplement `Jetson` as `Tegra`
- Bugfixes

## 1.1.0

### Changelog:

- Fix `Arduino`
- Extend Python interface

## 1.0.0

This is the initial release of Power Measurement Toolkit (PMT), as published in
[HUST'22](https://www.computer.org/csdl/proceedings-article/hust/2022/634900a044/1KnWyMBKdmo).

### Features:

- Comprehensive power measurement capabilities
- Support for Arduino and Python
- Initial toolkit structure and functionalities
