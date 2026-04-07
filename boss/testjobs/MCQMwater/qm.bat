@echo off
rem
echo.
echo ***********************************************************************
echo              MC Simulation with QM Solute
echo                 Single Solute in Water
echo                     25 C   1 atm
echo                  Equilibration 2.2M
echo                  Averaging     0-5M
echo             N = 500  T=25  P=1  Rcut = 10.0
echo ***********************************************************************
echo.
rem   For other solutes, replace qmzmat by the solute's Z-matrix,
rem   e.g., from the molecules/small directory.

set BOXES=%BOSSdir%
set BOSS=%BOSSdir%\boss
set LAMBDA=0.000 0.000 0.000

set INFILE=mcqmin
set UPFILE=mcqmup
set AVERAGE=mcqmav
set ZMATRIX=qmzmat
set SLVZMAT=slvzmat
set BANGPAR=%BOXES%\oplsaa.sb
set WATERBOX=%BOXES%\watbox
set ORG1BOX=%BOXES%\org1box
set ORG2BOX=%BOXES%\org2box
set SUMMARY=mcqmsum

rem     short NVT Equilibration

copy qmNVTpar + %BOXES%\oplsaa.par tmppar
SET PARAMETER=tmppar

set CONFIGURATIONS=200000
%BOSS% 111 %CONFIGURATIONS% %LAMBDA% -ot mcqmota -pl mcqmplta -sv mcqmsva

rem     Equilibration

del tmppar
copy qmpar + %BOXES%\oplsaa.par tmppar
set PARAMETER=tmppar

set CONFIGURATIONS=500000
%BOSS% 011 %CONFIGURATIONS% %LAMBDA% -ot mcqmota -pl mcqmplta -sv mcqmsva

for %%i in ( b c eq ) do %BOSS% 001 %CONFIGURATIONS% %LAMBDA% -ot mcqmot%%i -pl mcqmplt%%i -sv mcqmsv%%i

rem     Averaging

set CONFIGURATIONS=1000000
%BOSS% 011 %CONFIGURATIONS% %LAMBDA% -ot mcqmota -pl mcqmplta -sv mcqmsva

for %%i in ( b c d e5) do %BOSS% 001 %CONFIGURATIONS% %LAMBDA% -ot mcqmot%%i -pl mcqmplt%%i -sv mcqmsv%%i

del slvzmat
del tmppar
set BOXES=
set BOSS=
set CONFIGURATIONS=
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

