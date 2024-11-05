## Conversion of TCAD output for Pixelav simulation

The goal of this module is to take the output of TCAD simulation and prepare it for use with Pixelav simulation. The TCAD E field is used when simulating clusters.

The first step is to build the recipes_c-ansi library and compile the `gen_efield.c`
```
cd recipes_c-ansi
sudo make
cp include/nr*h /usr/local/include/ #For some reason, installation does not copy these headers
cd ..
./linkc_recipes gen_efield
```

The generation of the E field is done using a TCAD output folder, named `dot1_150x100_prd2022lk3_dj0305d` in this example. The two arguments are the folder name and the HV of the simulation for which to produce the field for. The script asks for two lines of input. The first are the number of gridpoints along each axis. The second line is the boundary where we force E field to be zero and is used to help the calculations converge. The boundary should be set to be larger than the silicon sensor thickness. Here we input 500 um.
```
./gen_efield dot1_150x100_prd2022lk3_dj0305d 475
25 17 92 #Input 1
500 #Input 2
```

The two outputs are `eplot.out` and `efield.out`. The former can be used for plotting by appending it to `ef_comp2024...` file to use with `tdr.py` in the `LA_trees_processing` module. The latter is used for simulation. Here we add four lines to the efield file:

```
Dot1_150x100_dj0305d@-475V,3.8T@90d,263K,rh=1.02/0.7,2.5/2.5,100Vdep
1.00 32600
3.80 0.0 0.00
285. 150. 100. 263. 2.5  2.5 1.02 0.70 0 1 25 17 92
```

* First line is the ascii header, used for bookkeeping
* Second line: momentum of particles to simulate. If it’s set to 1.0 it assumes 45 GeV, base run number for pixelav
* Third line: Magnetic field in tesla
* Fourth line: thickness of silicon, local x length, local y length, temperature in K, radiation exposures in units of neq 10^14 for electrons and holes (scale factors for trapping rates), hall factors for electrons and  holes, collected charge type (0=electrons, 1=holes), model for electron energy loss (we only use 1, it is from NIST), dimensions of E field array (we input it in previous step)
* Save it as pixel2.init file
* Run Pixelav simulation with new pixel2.init file

## Running TCAD
* Start from one of the folders that were used in previous TCAD simulations
* Simulation parameters are in files msh.grd, msh.dat and _des.cmd. _des.cmd is usually the one to modify, mostly the concentrations of two types of carrier traps.
* Run `dessis test_diode_dj0305d` #`test_diode_dj0305d` is the name of the folder containing details about TCAD simulation and it needs to be run from that folder
* Move the output folder to tcad_file_processing
