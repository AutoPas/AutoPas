# MD-Flexible

This programm aims to easily create and simulate molecular dynamic scenarios using AutoPas
## Compiling
build instructions for md-flexible
```bash
mkdir build
cd build
cmake [-G generator] [CMAKEOPTIONS] ..
ccmake [OPTIONS]..
make md-flexible [OPTIONS]
```
## Usage 
#####Information:
1. Options are eather inputed via a .yaml configuration file, via the command line or both 
2. Command line options overwrite the options specified in the yaml file
3. A sample YAML configuration file with all available options is located at: examples/md-flexible/input/FullConfigurationFile.yaml
5. For every execution, the options of the choosen scenario are printed on the command line 
#####Execution:
* Display all CL options
```bash
 examples/md-flexible/md-flexible --help
```
* Specify .yaml configuration file
```bash
 examples/md-flexible/md-flexible --yaml-filename <path-to-file>
```
* Default execution (3D Grid of 1000 particles)
```bash
 examples/md-flexible/md-flexible
```