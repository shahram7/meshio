@ECHO OFF
SETLOCAL enabledelayedexpansion

:: Set SIMULIA base directory here:
SET root=C:\SIMULIA\CAE

SET dir=%~dp0
SET plugindir=win_b64\code\python2.7\lib\abaqus_plugins
SET codedir=win_b64\code\bin

SET Index=1
FOR /d %%D IN (%root%\*) DO (
  SET "versions[!Index!]=%%D"
  SET /a Index+=1
)
set /a UBound=Index-2
for /l %%i in (1,1,%UBound%) do (
  ECHO "Installing Meshio-Plugin to !versions[%%i]!"
  SET str=!versions[%%i]!
  SET version=!str:~-4!
  :: Copy
  XCOPY %dir%meshio_plugin %root%\!version!\%plugindir%\MeshioPlugin\* /YQ
  XCOPY %dir%\..\..\meshio %root%\!version!\%codedir%\meshio\* /YQS

  XCOPY %dir%abq_meshio %root%\!version!\%codedir%\abq_meshio\* /YQ
)
PAUSE
