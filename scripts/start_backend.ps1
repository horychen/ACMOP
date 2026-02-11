# 使用 acmop.config.json 中配置的 conda 环境启动后端
# 默认环境名为 acmop
$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $PSScriptRoot
$ConfigPath = Join-Path $ProjectRoot "acmop.config.json"

$envName = "acmop"
if (Test-Path $ConfigPath) {
    try {
        $config = Get-Content $ConfigPath -Raw | ConvertFrom-Json
        $envName = $config.backend.virtualEnv.name
        if (-not $envName) { $envName = "acmop" }
    } catch {
        Write-Host "Warning: Could not parse acmop.config.json, using conda env: acmop"
    }
} else {
    Write-Host "No acmop.config.json found, using conda env: acmop"
}

$BackendDir = Join-Path $ProjectRoot "backend"
Set-Location $BackendDir
Write-Host "Starting backend with conda env: $envName"
& conda run -n $envName --no-capture-output uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
