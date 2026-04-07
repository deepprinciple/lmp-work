@echo off
rem
echo.
echo ***********************************************************************
echo            Potential Surface Scan for Anion with PDDG/PM3
echo ***********************************************************************
echo.
set BOXES=c:\boss
set BOSS=c:\boss\boss
rem the next line specifies the atom number for the scan            
set ATOM=000005
rem the first number below specifies the increment in Angstroms or degrees
rem between energy evaluations; next is 1., 2., or 3. for bond length,
rem angle or dihedral; then the total number of energy evaluations
set LAMBDA= 0.010 1.000 11.000

set INFILE=scanin
set UPFILE=scanup
set AVERAGE=scanav
set ZMATRIX=scanzmat
set SLVZMAT=scanzm
copy %BOXES%\scripts\pdgopt-par + %BOXES%\oplsaa.par scanpar
rem for Fletcher-Powell optimization, change
rem OPLSpar above to OPLSFpar
set PARAMETER=scanpar  
set BANGPAR=%BOXES%\oplsaa.sb
set WATERBOX=%BOXES%\watbox
set ORG1BOX=%BOXES%\org1box
set ORG2BOX=%BOXES%\org2box
set SUMMARY=scansum

rem  the 3 below requests the energy scan        
%BOSS% 311 %ATOM% %LAMBDA% -ot scanout -pl plt.pdb -sv scansva

del scanin
del scanpar
del scanzm
del scanup
del scanav
del scansva
set BOXES=
set BOSS=
set ATOM=
set LAMBDA=
set INFILE=
set UPFILE=
set AVERAGE=
set ZMATRIX=
set SLVZMAT=
set PARAMETER=
set BANGPAR=
set WATERBOX=
set ORG1BOX=
set ORG2BOX=
set SUMMARY=

