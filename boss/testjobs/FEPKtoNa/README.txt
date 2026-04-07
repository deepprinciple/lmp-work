
This is an FEP calculation for perturbing K+ to Na+ unbound
in methanol. As in test job fepcrown, 
the solvent consists of 255 OPLS-UA methanol molecules
from the stored box of 267.

The entire FEP is run by submission of just the f25.bat file.
Two windows are executed with delta lambda = +- 0.25.
Each window has 2M configurations of equilibration and 2M
of averaging. 10-A cutoffs are used. A larger system
size could be considered, though the
difference in free energy of binding would probably be affected
little as it is largely determined by relatively short-range
interactions, i.e., 1st and 2nd solvation shell.

The results from this job combine with those for the bound 
perturbation in test job fepcrown to yield the difference in
free energies of binding for the two ions in methanol.

Time required: 18 minutes on a 2.4 GHz Pentium.

Summary of results (easily obtained by executing xSUMDELG):
K+ -> Na+ Unbound
     lambda       DeltaG      Sigma       Sum
 0.0 to 0.25     -3.988       0.080      -3.988
 0.0 to 0.25     -4.314       0.083      -8.302
 0.0 to 0.25     -4.818       0.109     -13.120
 0.0 to 0.25     -5.027       0.085     -18.147

                      Sum of DeltaGs =  -18.147   
                      Sigma          =    0.180   

                 Exptl. DeltaG Solvation  -17.3

                   Computed DeltaG Bound  -16.415 +/-     0.153

            Computed DeltaDeltaG Binding    1.732 +/-     0.236
              Exptl. DeltaDeltaG Binding    2.45
                (JACS 104, 6895 (1982))


Other data from R. M. Izatt et al., Chem. Rev. 85, 271-339 (1985):

18-crown-6 at 25 C
                   water                    methanol
     ion     log K      delta G        log K      delta G

     Li+       0          0             0           0
     Na+      0.8       -1.09          4.36       -5.95
     K+       2.03      -2.77          6.08       -8.30
     Rb+      1.56      -2.13          5.32       -7.26
     Cs+      0.99      -1.35          4.71       -6.42
     NH4+     1.23      -1.68          4.14       -5.65
     CH3NH3+                           4.25       -5.80
