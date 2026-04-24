@echo off

if exist .env (
    for /f "usebackq tokens=* eol=#" %%i in (".env") do set %%i
)

if exist Scripts\python.exe (
    Scripts\python.exe -u pipexGUI.py 2>&1 | powershell -command "$input | Tee-Object -FilePath log.txt"
) else (
    python -u pipexGUI.py 2>&1 | powershell -command "$input | Tee-Object -FilePath log.txt"
)