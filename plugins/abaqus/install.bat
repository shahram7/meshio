@ECHO OFF
SETLOCAL enabledelayedexpansion

:: Set SIMULIA base directory here:
SET root=C:\SIMULIA\CAE\plugins

SET dir=%~dp0

SET Index=1
FOR /d %%D IN (%root%\*) DO (
  SET "versions[!Index!]=%%D"
  SET /a Index+=1
)
set /a UBound=Index-1
for /l %%i in (1,1,%UBound%) do (
  ECHO "Installing Meshio-Plugin to !versions[%%i]!"
  SET str=!versions[%%i]!
  SET version=!str:~-4!

  IF !version! GEQ 2020 (
    SET codedir=C:\SIMULIA\EstProducts\!version!\win_b64\code\bin
  ) ELSE (
    SET codedir=C:\SIMULIA\CAE\!version!\win_b64\code\bin
  )

  :: Copy
  XCOPY %dir%meshio_plugin %root%\!version!\MeshioPlugin\* /YQ

  XCOPY %dir%\..\..\meshio !codedir!\meshio\* /YQS
  XCOPY %dir%abq_meshio !codedir!\abq_meshio\* /YQ
)
PAUSE
