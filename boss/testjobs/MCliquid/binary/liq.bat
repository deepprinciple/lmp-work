@echo off
rem 
rem Command file for the BOSS program.  Execute by entering the following:
rem liq
echo.
echo ***********************************************************************
echo                 Binary Liquid Simulation
echo                   Equilibration 4000K 
echo                 Averaging     0000-4000K
echo                         T=25  P=1  
echo ***********************************************************************
echo.

set MCLIB=c:\boss
set boss=c:\boss\boss
set configurations=1000000
set lambda=0.0 0.0 0.0   

set INFILE=liqin
set UPFILE=liqup
set AVERAGE=liqav
set ZMATRIX=dummy.z
set SLVZMAT=liqzmat
copy liqpar + %MCLIB%\oplsaa.par tmppar
set PARAMETER=tmppar  
set BANGPAR=%MCLIB%\oplsaa.sb
set WATERBOX=%MCLIB%\watbox
set ORG1BOX=%MCLIB%\org1box
set ORG2BOX=%MCLIB%\org2box
set SUMMARY=liqsum

rem   Equilibration

rem   The first 1 below says start a new MC calculation from scratch.
rem   The solvent and solute coordinates are to be built.
rem   The second 1 says restart the averaging

%BOSS% 111 %CONFIGURATIONS% %LAMBDA% -ot liqota -pl liqplta -sv liqsva

for %%i in ( b c eq ) do %BOSS% 001 %CONFIGURATIONS% %LAMBDA% -ot liqot%%i -pl liqplt%%i -sv liqsv%%i
  
rem   Averaging

rem   The 01 below says continue the MC run, but restart the averaging.

%BOSS% 011 %CONFIGURATIONS% %LAMBDA% -ot liqotd -pl liqpltd -sv liqsvd

for %%i in ( e f e4 ) do %BOSS% 001 %CONFIGURATIONS% %LAMBDA% -ot liqot%%i -pl liqplt%%i -sv liqsv%%i

del tmppar liq*b liq*c liq*d liq*e liq*f
set MCLIB=
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

