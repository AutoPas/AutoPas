# MD-Flexible 

This demo shows how to easily create and simulate molecular dynamic
scenarios using AutoPas and flexibly configure all of it's options.

## Compiling
To build MD-Flexible execute the following from the AutoPas root folder:
```bash
mkdir build && $_
cmake ..
make md-flexible
```

## Usage 

When running md-flexible without any arguments a default simulation with
all AutoPas options active is run and it's configuration printed. From
there you can restrict any AutoPas options or change the simulation.

For all available option see:
```bash
 examples/md-flexible/md-flexible --help
```

### Input

MD-Flexible accepts input via command line arguments and YAML files.
When given both, any command line options will overwrite their YAML
counterparts.

### Checkpoints

MD-Flexible can be initialized through a previously written VTK file.
Please use only VTK files written by MD-Flexible since the parsing is
rather simple. The VTK file only contains Information about all
particles positions, velocities, forces and typeIDs. All other options,
especially the simulation box size and particle properties (still) need
to be set through a YAML file.
