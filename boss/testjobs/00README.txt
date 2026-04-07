Most directories contain a file named README or 00README
that should be read.


                BOSS Test Jobs

To learn basic operations of BOSS, it is recommended
that the test jobs be reviewed and processed in the
following order. In all cases, it is important to
look carefully at the Z-matrix, par and cmd (or bat) files.
The point is to understand the test jobs, not just run
them. Descriptions of the test jobs are in the user's manual.

Basic Molecular & Quantum Mechanics
1 - optimize
2 - consearch
3 - dihdrive
4 - scan

Standard MC Simulations
5 - MCgas  - molecule in the gas phase
6 - MCGBSA - molecule with GB/SA hydration
7 - MCwater - liquid water 
8 - MCdroplet - water droplet
9 - MCliquid - various pure liquids 
10- ionwater - an ion in water
11- MCsolute - long MC runs for a solute in TIP4P water
12- CustomSolvent - examples for a solute in a custom solvent.

MC/FEP Simulations
13- FEPgas (many basic X to Y examples in the gas phase)
14- FEPaq  (many basic X to Y examples in a water box)
15- SOSgas - uses overlap sampling (")
16- SOSaq  - uses overlap sampling (")
17- FEPKtoNa - FEP for K+ to Na+ in methanol
18- FEPcrown - " " "  "  bound to 18-crown-6
19- FEPannihil - FEP to compute free energies of hydration
20- dihpmf - FEP for a torsion angle
21- pairpmf - FEP for separating two species
22- rxnpmf - FEP for a reaction profile - QM/MM
23- FEPpolrz - calculations with solvent polarization

More Specialized Simulations
24- MCslab
25- ionnonaq
26- MUSIC - like MCSS - saturate protein with probe molecules
27- FEPcap
28- FEPflex
29- linres - organic solute in water & estimate properties
30- MCQMwater - MC simulation for a QM solute in TIP4P water
31- HalogenBond - treatment of halogens with X sites
32- cationpi - optimizations for cation-pi complexes inc. 1/r4 terms
33- MixSolvent - octanol/water mixed solvent
