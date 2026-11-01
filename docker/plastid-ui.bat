@echo off
rem PlastidToolSuite on Windows. Double-click this file. Your work goes in Documents\PlastidToolSuite.
setlocal
set IMAGE=plastid-suite
set WORKSPACE=%USERPROFILE%\Documents\PlastidToolSuite
if not exist "%WORKSPACE%" mkdir "%WORKSPACE%"

docker version >nul 2>&1
if errorlevel 1 (
  echo Docker Desktop is not running. Start it, wait until it says it is running, then try again.
  pause
  exit /b 1
)

docker image inspect %IMAGE% >nul 2>&1
if errorlevel 1 (
  echo Building the toolkit. This happens once, needs the network, and takes a long while.
  docker build -f "%~dp0Dockerfile" -t %IMAGE% "%~dp0.."
  if errorlevel 1 (
    echo The build stopped. The lines above say why.
    pause
    exit /b 1
  )
)

for /f %%t in ('powershell -NoProfile -Command "[guid]::NewGuid().ToString('N')"') do set TOKEN=%%t
start /b cmd /c "timeout /t 5 >nul & start http://127.0.0.1:26214/?token=%TOKEN%"
echo Keep this window open while you work. Close it to stop the interface.
docker run --rm -p 127.0.0.1:26214:26214 -e PLASTID_UI_TOKEN=%TOKEN% -v "%WORKSPACE%:/work" %IMAGE%
