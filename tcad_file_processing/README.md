##Conversion of TCAD output for Pixelav simulation

The goal of this module is to take the output of TCAD simulation and prepare it for use with Pixelav simulation. The TCAD E field is used when simulating clusters.

The first step is to build the recipes_c-ansi library and compile the `gen_efield.c`
```
cd recipes_c-ansi
sudo make
cp include/nr*h /usr/local/include/ #For some reason, installation does not copy these headers
cd ..
./linkc_recipes gen_efield
```

The generation of the E field is done using a TCAD output folder, named `dot1_150x100_prd2022lk3_dj0305d` in this example. The two arguments are the folder name and the HV of the simulation for which to produce the field for. The script asks for two lines of input. The first are the number of gridpoints along each axis. The second line is something that should be larger than the silicon sensor thickness (I have to ask Morris what this is exactly), we input 500 um.
```
./gen_efield dot1_150x100_prd2022lk3_dj0305d 475
25 17 92 #Input 1
500 #Input 2
```