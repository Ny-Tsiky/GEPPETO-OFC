@echo off

call "setup_mingw.bat"

cd .

chcp 1252

if "%1"=="" ("D:\matlab2021\bin\win64\gmake"  -f my_predict_rtw.mk all) else ("D:\matlab2021\bin\win64\gmake"  -f my_predict_rtw.mk %1)
@if errorlevel 1 goto error_exit

exit /B 0

:error_exit
echo The make command returned an error of %errorlevel%
exit /B 1