@echo off
rem 
rem Command file for the BOSS program.  Execute by entering the following:
rem gas
echo.
echo ***********************************************************************
echo                 Gas-phase Monte Carlo
echo                   Equilibration 1.0 M 
echo                       Averaging 2.0 M 
echo ***********************************************************************
echo.

set boss=%BOSSdir%\boss
set configurations=200000
set lambda=0.0 0.0 0.0   

set INFILE=gasin
set UPFILE=gasup
set AVERAGE=gasav
set ZMATRIX=gaszmat
set SLVZMAT=slvzmat
copy gaspar + %BOSSdir%\oplsaa.par tmppar
set PARAMETER=tmppar  
set BANGPAR=%BOSSdir%\oplsaa.sb
set WATERBOX=%BOSSdir%\watbox
set ORG1BOX=%BOSSdir%\org1box
set ORG2BOX=%BOSSdir%\org2box
set SUMMARY=gassum

rem  Equilibration

%BOSS% 111 %CONFIGURATIONS% %LAMBDA% -ot gasota -pl gasplta -sv gassva

for %%i in ( b c d e ) do %BOSS% 001 %CONFIGURATIONS% %LAMBDA% -ot gasot%%i -pl gasplt%%i -sv gassv%%i

rem   Averaging 

set configurations=400000
%BOSS% 011 %CONFIGURATIONS% %LAMBDA% -ot gasota -pl gasplta -sv gassva

for %%i in ( b c d e ) do %BOSS% 001 %CONFIGURATIONS% %LAMBDA% -ot gasot%%i -pl gasplt%%i -sv gassv%%i

del slvzmat
del gassv*
if exist gasout del gasout
ren gasote gasout
del gasot*
ren gasplte gplt.pdb
del gasplt*
ren gplt.pdb gasplt.pdb
del tmppar
set boss=
set configurations=
set lambda=
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

