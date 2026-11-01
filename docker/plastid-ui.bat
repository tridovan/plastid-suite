@echo off
rem PlastidToolSuite on Windows. Double-click this file. Your work goes in Documents\PlastidToolSuite.
setlocal
set IMAGE=plastid-suite
set NAME=plastid-ui
set PORT=26214
set WORKSPACE=%USERPROFILE%\Documents\PlastidToolSuite
if not exist "%WORKSPACE%" mkdir "%WORKSPACE%"

docker version >nul 2>&1
if errorlevel 1 (
  echo Docker Desktop is not running. Start it, wait until it says it is running, then try again.
  pause
  exit /b 1
)

rem A window closed earlier can leave the toolkit running. Open that one rather than start a second.
for /f %%c in ('docker ps -q --filter "name=^%NAME%$"') do goto :reopen

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
rem The browser opens once the server answers: 403 without the session cookie is an answer.
start /b powershell -NoProfile -Command "for($i=0;$i -lt 120;$i++){try{Invoke-WebRequest -UseBasicParsing -TimeoutSec 2 'http://127.0.0.1:%PORT%/' | Out-Null; break}catch{if($_.Exception.Response){break}; Start-Sleep 1}}; Start-Process 'http://127.0.0.1:%PORT%/?token=%TOKEN%'"
docker run --rm --init --name %NAME% -p 127.0.0.1:%PORT%:%PORT% -e PLASTID_UI_TOKEN=%TOKEN% -v "%WORKSPACE%:/work" %IMAGE%
if errorlevel 1 (
  echo The toolkit stopped. If that was not you, the lines above say why.
  pause
)
exit /b 0

:reopen
for /f %%t in ('docker exec %NAME% printenv PLASTID_UI_TOKEN') do set TOKEN=%%t
start "" "http://127.0.0.1:%PORT%/?token=%TOKEN%"
exit /b 0
