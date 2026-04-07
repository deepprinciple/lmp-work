@echo off
rem 
rem Command file for the BOSS program.  Execute by entering the following:
rem liq
echo.
echo ***********************************************************************
echo                 Pure Liquid Simulation
echo                 Restart Averaging and Run 4M
echo ***********************************************************************
echo.

set boss=%BOSSdir%\boss
set configurations=500000
set lambda=0.0 0.0 0.0   

set INFILE=liqin
set UPFILE=liqup
set AVERAGE=liqav
set ZMATRIX=dummy.z
set SLVZMAT=liqzmat
copy liqpar + %BOSSdir%\oplsaa.par tmppar
set PARAMETER=tmppar  
set BANGPAR=%BOSSdir%\oplsaa.sb
set WATERBOX=%BOSSdir%\watbox
set ORG1BOX=%BOSSdir%\org1box
set ORG2BOX=%BOSSdir%\org2box
set SUMMARY=liqsum

rem   Averaging

rem   The 01 below says continue the MC run, but restart the averaging.

%BOSS% 011 %CONFIGURATIONS% %LAMBDA% -ot liqota -pl liqplta -sv liqsva

for %%i in ( b c e2 a b c e4 ) do %BOSS% 001 %CONFIGURATIONS% %LAMBDA% -ot liqot%%i -pl liqplt%%i -sv liqsv%%i

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

