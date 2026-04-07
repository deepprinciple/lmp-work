@echo off
rem     usage: xPDGOPT- molecule
rem     file molecule.z   must exist

echo    *****************************
echo     PDDG/PM3 Anion Optimization
echo              %1
echo    *****************************

copy %1.z optzmat
rem   optimize geometry
call c:\boss\scripts\PDGOPT-

echo
echo    Finished

