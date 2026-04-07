@echo off
rem
rem Script to run MC calculations at multiple pressures
rem
rem      usage:      xLIQP 
rem
call liq
cd P2K
copy ..\dummy.z .
call liq2K
cd ..
cd P5K
copy ..\dummy.z .
call liq5K
cd ..
cd P10K
copy ..\dummy.z .
call liq10K
