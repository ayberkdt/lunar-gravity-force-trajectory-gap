param(
    [Parameter(Mandatory=$true)]
    [string]$Deadline
)

$ErrorActionPreference = 'Continue'
$pythonExe = 'D:\Masaustu\LUNAR_SIMULATION\.venv\Scripts\python.exe'
$codeDir = 'C:\Users\ayber\Desktop\Makale\codebase\python_codes'
$metricsDir = 'C:\Users\ayber\Desktop\Makale\codebase\metrics'
$blendScript = Join-Path $codeDir 'rev10_blend_lro_convergence.py'
$convergenceScript = Join-Path $codeDir 'rev10_sobol_convergence.py'
$manifestPath = Join-Path $metricsDir 'r10_overnight_queue_manifest.json'
$deadlineTime = [DateTimeOffset]::Parse($Deadline)
$started = [DateTimeOffset]::Now
$steps = @()

function Add-StepRecord {
    param([string]$Name, [DateTimeOffset]$Start, [DateTimeOffset]$End, [int]$ExitCode)
    $script:steps += [pscustomobject]@{
        name = $Name
        started_utc = $Start.ToUniversalTime().ToString('o')
        ended_utc = $End.ToUniversalTime().ToString('o')
        wall_s = ($End - $Start).TotalSeconds
        exit_code = $ExitCode
    }
}

Write-Output "[queue] started=$($started.ToString('o')) deadline=$($deadlineTime.ToString('o'))"

if ([DateTimeOffset]::Now -lt $deadlineTime) {
    $stepStart = [DateTimeOffset]::Now
    Write-Output '[queue] step 1: 28-day LRO corrected-blend convergence'
    & $pythonExe -u $blendScript run --deadline $Deadline
    $stepCode = $LASTEXITCODE
    $stepEnd = [DateTimeOffset]::Now
    Add-StepRecord -Name 'blend_lro_28day_convergence' -Start $stepStart -End $stepEnd -ExitCode $stepCode
    Write-Output "[queue] step 1 exit=$stepCode wall_s=$(($stepEnd-$stepStart).TotalSeconds)"
}

if ([DateTimeOffset]::Now -lt $deadlineTime) {
    $stepStart = [DateTimeOffset]::Now
    Write-Output '[queue] step 2: 17-orbit selective Sobol convergence'
    & $pythonExe -u $convergenceScript run --deadline $Deadline
    $stepCode = $LASTEXITCODE
    $stepEnd = [DateTimeOffset]::Now
    Add-StepRecord -Name 'sobol_selective_convergence' -Start $stepStart -End $stepEnd -ExitCode $stepCode
    Write-Output "[queue] step 2 exit=$stepCode wall_s=$(($stepEnd-$stepStart).TotalSeconds)"
}

$blendStatus = $null
$convergenceStatus = $null
$blendPath = Join-Path $metricsDir 'r10_blend_lro_convergence.json'
$convergencePath = Join-Path $metricsDir 'r10_sobolA_convergence.json'
if (Test-Path $blendPath) {
    $blend = Get-Content -Raw $blendPath | ConvertFrom-Json
    $blendStatus = [pscustomobject]@{
        complete = $blend.complete
        record_count = $blend.records.Count
        skipped_for_deadline = $blend.skipped_for_deadline
        timing_comparable = $blend.timing_comparable
        summary_available = $null -ne $blend.summary
    }
}
if (Test-Path $convergencePath) {
    $convergence = Get-Content -Raw $convergencePath | ConvertFrom-Json
    $convergenceStatus = [pscustomobject]@{
        complete = $convergence.complete
        completed_orbits = $convergence.rows.Count
        selected_count = $convergence.selected_count
        stopped_for_deadline = $convergence.stopped_for_deadline
        timing_comparable = $convergence.timing_comparable
        summary = $convergence.summary
    }
}

$ended = [DateTimeOffset]::Now
$manifest = [ordered]@{
    schema = 'r10_overnight_queue_manifest_v1'
    started_utc = $started.ToUniversalTime().ToString('o')
    deadline = $deadlineTime.ToString('o')
    ended_utc = $ended.ToUniversalTime().ToString('o')
    total_wall_s = ($ended - $started).TotalSeconds
    steps = $steps
    blend_status = $blendStatus
    convergence_status = $convergenceStatus
}
$json = $manifest | ConvertTo-Json -Depth 10
[System.IO.File]::WriteAllText($manifestPath, $json, [System.Text.UTF8Encoding]::new($false))
Write-Output "[queue] manifest=$manifestPath"
Write-Output "[queue] finished=$($ended.ToString('o'))"
