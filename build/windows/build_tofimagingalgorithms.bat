call set_vars.bat
pushd .

set REPOSPATH=%WORKSPACE%\ToFImaging
set DEST=%WORKSPACE%\builds
set SPECSTR=-spec win32-msvc

mkdir %DEST%\build-tofimagingalgorithm
cd %DEST%\build-tofimagingalgorithm

%QTBINPATH%\qmake.exe -makefile -r %SPECSTR% ..\..\ToFImaging\TOF_ImagingAlgorithm\TOF_ImagingAlgorithm.pro -o Makefile
%QTBINPATH%\..\..\..\Tools\QtCreator\bin\jom.exe -f Makefile clean
%QTBINPATH%\..\..\..\Tools\QtCreator\bin\jom.exe -f Makefile mocables all
%QTBINPATH%\..\..\..\Tools\QtCreator\bin\jom.exe -f Makefile release

echo "Build tests"

cd %REPOSPATH%\UnitTests\

FOR /D %%I IN (*) DO @call :testloopbody %REPOSPATH% %%I %DEST%

echo "Tests built"
popd 

goto :eof

:testloopbody
echo %2
if exist "%1\UnitTests\%2\%2.pro" (
	mkdir %3\%2
	cd %3\%2

	%QTBINPATH%\qmake.exe -makefile ..\..\ToFImaging\UnitTests\%2\%2.pro -o Makefile
	%QTBINPATH%\..\..\..\Tools\QtCreator\bin\jom.exe -f Makefile clean
	%QTBINPATH%\..\..\..\Tools\QtCreator\bin\jom.exe -f Makefile mocables all
	%QTBINPATH%\..\..\..\Tools\QtCreator\bin\jom.exe -f Makefile release
)
goto :eof