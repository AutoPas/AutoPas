# MD-Flexible

This programm aims to easily create and simulate different molecular dynamic scenarios using AutoPas
## Compiling
1. Create build directory: <br />
mkdir build<br />
cd build<br />
2. Run CMake:<br />
cmake [-G generator] [CMAKEOPTIONS] ..<br />
3. Run the ccmake command to select further build options:<br />
ccmake [OPTIONS]..<br />
3. Building:<br />
make md-flexible [OPTIONS]
## Usage
For every execution, the options of the choosen scenario are printed on the command line <br />
Options are eather inputed via a .yaml configuration file, via the command line or both <br />
Command line options overwrite the options specified in the yaml file <br />
Use --help to display all CL options <br />
To specify the .yaml file, use --yaml-filename <br />
A sample .yaml configuration file with all available options is located at: <br /> examples/md-flexible/input/FullConfigurationFile.yaml

####Execution: 
 examples/md-flexible/md-flexible [OPTIONS]
