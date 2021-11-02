# MD-Flexible

This demo shows how to easily create and simulate molecular dynamic
scenarios using AutoPas and flexibly configure all of it's options.

## Documentation
The documentation can be found at our website:
 <https://autopas.github.io/doxygen_documentation_md-flexible/git-master/index.html>

Alternatively you can build the documentation on your own:
* requirements: [Doxygen](http://www.doxygen.nl/)
* `make doc_doxygen_md-flexible`

## Compiling
To build MD-Flexible execute the following from the AutoPas root folder:
```bash
mkdir build && cd $_
cmake ..
make md-flexible
```

### Compiling with MPI
To use the MPI parallelization of md-flexible activate `MD_FLEXIBLE_USE_MPI` via `cmake`:
```bash
cmake -DMD_FLEXIBLE_USE_MPI=ON ..
Using this option does not guarantee that MPI will be used. There are some additional requirements:
* MPI is installed
* mpirun is called with more than 1 process
```

## Testing
Simple tests can be run via:
```bash
make mdFlexTests
ctest -R mdFlexTests
```
To run Tests with MPI enabled:
```bash
mpiexec -np 4 ./examples/md-flexible/tests/mdFlexTests
```

## Usage

When starting md-flexible without any arguments a default simulation with
most AutoPas options active is run and it's configuration printed. From
there you can change any AutoPas options or change the simulation.

For all available arguments and options see:
```bash
 examples/md-flexible/md-flexible --help
```

The Falling Drop example simulation can be started with:
```bash
 cd examples/md-flexible
 ./md-flexible --yaml-filename fallingDrop.yaml
```

To execute md-flexible using MPI run the following command in the build directory:
```bash
mpirun -np 4 ./examples/md-flexible/md-flexible --yaml-filename ./examples/md-flexible/fallingDrop.yaml
```
The number 4 can be exchanged by any number assuming the hardware supports the 
selected number of processes.

### Input

MD-Flexible accepts input via command line arguments and YAML files.
When given both, any command line options will overwrite their YAML
counterparts.

The keywords for the YAML file are the same as for the command line
input. However, since there are options that can only be defined
through the YAML file there is also the file [`input/AllOptions.yaml`](https://github.com/AutoPas/AutoPas/blob/master/examples/md-flexible/input/AllOptions.yaml)
to be used as a reference.

Additionally, that options that require a list of choices also
accept the keyword `all` to enable all, even discouraged choices.

#### Object generators

To quickly set up scenarios md-flexible provides a couple of object
generators that create 3D shapes filled with particles. From the command line
only one generator can be used at a time, however when using a YAML file one
can use as an arbitrary amount of generators. In YAML files it is also
possible to generate multiple objects from the same generator as one
specifies the objects directly. For a list of all possible objects and their
descriptions see [`src/Objects`](https://autopas.github.io/doxygen_documentation_md-flexible/git-master/dir_8e5023335c6d80afeb9fe41ac1daf95f.html).
For examples how to define and configure each object see [`input/AllOptions.yaml`](https://github.com/AutoPas/AutoPas/blob/master/examples/md-flexible/input/AllOptions.yaml).

### Output

* After every execution, a configuration YAML file is generated. It is possible
  to use this file as input for a new simulation.

* You can generate vtk output by providing a vtk-filename
(see help for details). The output contains two different vtk files. One for
  visualizing the particles another for the visualization of the subdomains of
  the domain decomposition.
The cells contain additional information about the configuration of the AutoPas
  container responsible for simulating the cell.
To visualize the particle records load the .pvtu files in ParaView. To visualize
  the cells of the decomposition load the .pvts files in ParaView.


### Checkpoints

MD-Flexible can be initialized through a previously written VTK file.
Please use only VTK files written by MD-Flexible since the parsing is
rather strict. The VTK file only contains Information about all
particles positions, velocities, forces and typeIDs. All other options,
especially the simulation box size and particle properties (still) need
to be set through a YAML file.

### Command line Completions

md-flexible can generate a shell completions file with it's latest options.
Feel free to add completions for your favorite shell.

#### zsh

1. Run:
```zsh
./md-flexible --zsh-completions
```
This will generate a file `_md-flexible` containing the completions definitions. 
Store it where you feel appropriate.
 
2. In your `.zshrc` prepend the path to the directory where `_md-flexible` is saved to `fpath`:
```zsh
fpath=("HereGoesThePath" $fpath)
```

3. Initialize the auto complete system by adding this line to `.zshrc` after the `fpath` line.:
```zsh
autoload -U compinit && compinit
```

4. Reload your `zsh`

### Misc

* `tuning-phases` overwrites restrictions set by `iterations`
