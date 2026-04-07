@echo off
rem 
rem Command file for the BOSS program.  Execute by entering the following:
rem liq
echo.
echo ***********************************************************************
echo              Flexible Pure Liquid Simulation
echo                     Equilibration 5M    
echo                       Averaging 10M          
echo                        T=25  P=10K                 
echo ***********************************************************************
echo.

set BOSS=%BOSSdir%\boss
set configurations=1000000
set lambda=0.0 0.0 0.0   

copy ..\P5K\liqin .
set INFILE=liqin
set UPFILE=liqup
set AVERAGE=liqav
set ZMATRIX=dummy.z
set SLVZMAT=liqzmat
copy liqpar10K + %BOSSDIR%\oplsaa.par tmppar
set PARAMETER=tmppar  
set BANGPAR=%BOSSDIR%\oplsaa.sb
set WATERBOX=%BOSSDIR%\watbox
set ORG1BOX=%BOSSDIR%\org1box
set ORG2BOX=%BOSSDIR%\org2box
set SUMMARY=liqsum

rem   Equilibration

rem   The first 1 below says start a new MC calculation from scratch.
rem   The solvent and solute coordinates are to be built.
rem   The second 1 says restart the averaging

%BOSS% 011 %CONFIGURATIONS% %LAMBDA% -ot liqota -pl liqplta -sv liqsva

for %%i in ( b c d eq ) do %BOSS% 001 %CONFIGURATIONS% %LAMBDA% -ot liqot%%i -pl liqplt%%i -sv liqsv%%i
  
rem   Averaging

rem   The 01 below says continue the MC run, but restart the averaging.

%BOSS% 011 %CONFIGURATIONS% %LAMBDA% -ot liqota -pl liqplta -sv liqsva

for %%i in ( b c d e5 a b c d e10 ) do %BOSS% 001 %CONFIGURATIONS% %LAMBDA% -ot liqot%%i -pl liqplt%%i -sv liqsv%%i

del tmppar
set BOSS=
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

