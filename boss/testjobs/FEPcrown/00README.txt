
This is an FEP calculation for perturbing K+ to Na+ bound
to 18-crown-6 in methanol. OPLS-AA is used for the solutes
and the solvent consists of 255 OPLS-UA methanol molecules
from the stored box of 267. 10-A cutoffs are used.

The entire FEP is run by submission of just the f25.bat file.
Two windows are executed with delta lambda = +- 0.25.
Each window has 2M configurations of equilibration and 2M
of averaging. The calculations could be extended and more
windows could be run for more precise results.

The results from this job combine with those for the unbound 
perturbation in test job fepktona to yield the difference in
free energies of binding for the two ions in methanol.

K+ -> Na+ bound

Time required:  28  min. on a  2.4 GHz Pentium IV.

Summary of results (obtained easily with xSUMDELG):
     lambda       DeltaG,      Sigma,      Sum
 0.0 -> 0.25      -4.435,      0.058,     -4.435
0.25 -> 0.50      -4.007,      0.037,     -8.442
0.50 -> 0.75      -4.028,      0.121,    -12.470
0.75 -> 1.00      -3.945,      0.064,    -16.415

 Sum of DeltaGs =  -16.414999
 Sigma          =    0.153199

See the README file in fepktona for the net binding results.

valala contains an FEP calculation for converting valine to alanine bound to 18-c-6.
