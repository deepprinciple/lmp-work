
MC simulations of liquid TIP4P water with 512 water molecules. 
The calculations start from the stored, cubic water box.
T = 25 C, P = 1 atm, and the OO cutoff distance is 9 A in all cases.

cube contains the setup and results for a cubic periodic cell (normal).

rhombic is for an orthorhombic periodic cell.

slab is for a water slab, periodic in x and y with vacuum on the +z and -z faces.

The files in the three directories are the same except for designation of
SLAB and VXYZ on line 9 in the t4ppar files.

t4p.bat or t4pcmd is the command file for the initial run with 5M 
configurations of equibration and 5M of averaging.
t4prest.bat shows how to restart from there and run 20M configurations of averaging.
t4pcont.bat continues for an additional 20M configurations of averaging to yield
a total of 40M.
Each run of 2M configurations takes about 7 minutes on a 2.4 GHz Pentium IV.

Prior output files at the end of the 40M configurations of averaging are called ote40.

Note: One can use the same files to run simulations for the
other stored liquids by simply changing line 3 of the t4ppar file
to the name of the other liquid (the names are listed under
SVMOD1 in section 9 of the user's manual). Also change the number
of solvent molecules from 512 to 267.
