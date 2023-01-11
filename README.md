# glyme-peo-analysis

## Overview

Program developed for analysis of Molecular Dynamics simulations involving systems with salts dissolved in solvents able to form coordinational bonds. 
For cations coordination numbers with details about number of solvent/anion atoms and solvent/anion molecules engaged in coordination are calculated.
In addition it is possible to calculate residence time autocorrelation functions and Venn diagrams for coordination bonds and determine for anions/solvent molecules number of coordinated cations for each time step. 

For reference see: https://doi.org/10.1021/acs.jpcb.1c05793

## Requirements and compilation

Only C language compiler (here GCC) and GNU Make are required.

To compile simply run
```
make
```

## Running

For the simplest case with default 
```
./program.x [ options ] xyz_file system_info_file
```

As an input the program needs trajectory from MD in xyz format and file with description of molecules present in the system (described below).

Options:
* `-h` - print help
* `-s value` - solvent threshold - threshold for determination in Angstrom of coordination between cation and solvent atom, default 1.5
* `-a value` - anion threshold - analogously as above, threshold in Angstrom for coordination between cation and anion, default 1.5
* `-b value` - simulation box size in Angstrom, default 22.00, note that only cubic systems are considered so only one value is given
* `-o` - produce one output file for all cations instead of separate file for every cation
* `-r` - calculate residence functions for solvent molecules
* `-f` - calculate residence functions for anions
* `-v` - calculate data for Venn diagrams
* `-d` - create files with coordination data for solvent molecules
* `-p name` - file name with box sizes for each timestep in NpT ensemble (without this option it is assumed that simulations was in NVT ensemble)

## System info file

As an input the program needs a file with description of molecules in xyz file in the following format:
* first line containing number of defined entries N
* N lines in the format:
    * type of entry: `cation/anion/solvent`
    * name of entry (max 20 characters)
    * the symbol of the first atom of the entry (max 4 characters)
    * number of molecules of this type
    * number of atoms in the molecule
    * symbol of atom considered as the one side of coordination bond
    * number of such atoms in one molecule
 
Note that entries should be listed in the same order as they appear in the xyz file.
 
For example, in the system containing 30 ionic pairs Na-TFSI and 300 monoglyme molecules in which we want to look for the Na - O(TFSI/monoglyme) bonds, where in the xyz file at the beginning are listed positions of Na cations, then TFSI anions and in the end solvent molecules, this file has the following form:
```
3
cation Na Na 30 1 Na 1
anion TFSI O 30 30 O 4
solvent monoglyme C 300 22 O 2
```

This file with example xyz file are located in `samples` directory.

## Output description

For provided samples with running:
```
./program.x -s 2.5 -a 2.5 -b 39.1619 -r -f -v -d samples/m-t.39.1619.xyz samples/m-t-sys-info.txt
```

There are produced some output files:
* `Na_x.dat` - files containing in every line data in format: 
```
step_number coordination_number corrdinating_atoms_from_solvent coordinating_solvent_molecules coordinating_atoms_from_anion coordinating_anions
```
* `monoglyme_x.dat` - output file for each solvent molecule with every line in format:
```
step_number number_of_distinct_coordinated_cations number_of_solvent_atoms_involved_in_coordination
```
* `cations_anions.dat` - file containing in every line number of currently coordinated cations to each anion, data for creating histogram
* `anions_times.dat` - file containing number of frames with existence of coordination between anion and cation, data for creating histogram
* `solvent-times.dat` - as above but for solvent molecules
* `residence-times-TFSI.dat` - residence autocorrelation function for bonding with anions without distinguishing particular anion atoms, i.e. when O1 atom from given TFSI anion forms at time `t` bond with given cation and in time `t + 1` no, but O2 atom forms a bond with the same cation it is treaten as non-broken between these two timesteps
* `residence-times-atoms-TFSI.dat` - as above, but treating each tracked atom of anion separately, i.e. when O1 atom from given TFSI anion forms at time `t` bond with given cation and in time `t + 1` no, but O2 atom forms a bond with the same cation it is treaten as broken between these two timesteps
* `residence-times-monoglyme.dat`, `residence-times-atoms-monoglyme.dat` - as above but for solvent (monoglyme)
* `venn-cations-monoglyme.dat` - data for Venn diagram from cations perspective
* `venn-solvent-monoglyme.dat` - as above, but from solvent molecules perspective
