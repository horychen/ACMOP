# Start Frontend and Backend for ACMOP
$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $PSScriptRoot
$ConfigPath = Join-Path $ProjectRoot "acmop.config.json"

# Default configuration
$envName = "acmop"
$BackendFolder = "backend_v2"

# Try to parse configuration
if (Test-Path $ConfigPath) {
    try {
        $config = Get-Content $ConfigPath -Raw | ConvertFrom-Json
        if ($config.backend.virtualEnv.name) {
            $envName = $config.backend.virtualEnv.name
        }
    } catch {
        Write-Host "Warning: Could not parse acmop.config.json, using conda env: $envName"
    }
} else {
    Write-Host "No acmop.config.json found, using conda env: $envName"
}

# Define Paths
$FrontendDir = Join-Path $ProjectRoot "frontend"

# NOTE: Modify "backend" to "backend_v2" if you are currently using backend_v2
$BackendDir = Join-Path $ProjectRoot $BackendFolder

Write-Host "Starting Frontend in a new window..."
# Starts frontend in a new PowerShell window
Start-Process powershell -ArgumentList "-NoExit", "-Command", "Set-Location '$FrontendDir'; npm run dev"

Write-Host "Starting Backend ($BackendFolder) in a new window using conda env: $envName..."
# Starts backend in a new PowerShell window
$BackendCmd = "Set-Location '$BackendDir'; conda run -n $envName --no-capture-output uvicorn main:app --reload --host 0.0.0.0 --port 8000"
Start-Process powershell -ArgumentList "-NoExit", "-Command", $BackendCmd

Write-Host "Both services have been launched in separate windows."
