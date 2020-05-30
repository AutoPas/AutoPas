# MD-Flexible

This demo shows how to easily create and simulate molecular dynamic
scenarios using AutoPas and flexibly configure all of it's options.

## Documentation
The documentation can be found at our website:
 <https://www5.in.tum.de/AutoPas/doc_doxygen_md-flexible/master/>

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

## Usage

When starting md-flexible without any arguments a default simulation with
most AutoPas options active is run and it's configuration printed. From
there you can change any AutoPas options or change the simulation.

For all available arguments and options see:
```bash
 examples/md-flexible/md-flexible --help
```

Have a look at the `completion/` directory for shell auto completion scripts.

### Input

MD-Flexible accepts input via command line arguments and YAML files.
When given both, any command line options will overwrite their YAML
counterparts.

The keywords for the YAML file are the same as for the command line
input. However, since there are options that can only be defined
through the YAML file there is also the file `input/ALLOptions.yaml`
to be used as a reference.

### Output

* After every execution, a configuration YAML file is generated. It is
possible to use this file as input for a new simulation.
* You can generate vtk output by providing a vtk-filename
(see help for details).

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
