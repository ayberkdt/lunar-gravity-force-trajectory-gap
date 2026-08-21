# R44 overnight runner: all eight registered cells, sequentially, 10 workers.
# Launch ONLY on the user's explicit go signal. Estimated ~3 h total.
# Usage:  powershell -ExecutionPolicy Bypass -File run_r44_overnight.ps1
$ErrorActionPreference = "Continue"
Set-Location $PSScriptRoot
$cells = @(
    @{d = "A"; b = "1.00"}, @{d = "B"; b = "1.00"},
    @{d = "A"; b = "0.75"}, @{d = "B"; b = "0.75"},
    @{d = "A"; b = "0.50"}, @{d = "B"; b = "0.50"},
    @{d = "A"; b = "1.25"}, @{d = "B"; b = "1.50"}
)
$t0 = Get-Date
foreach ($c in $cells) {
    $tag = "$($c.d)@$($c.b)"
    Write-Output "=== R44 $tag start $(Get-Date -Format HH:mm:ss) ==="
    python rev44_equal_work_tighter.py run --design $c.d --beta $c.b `
        --workers 10 --deadline-min 60 2>&1 |
        Tee-Object -Append -FilePath ..\metrics\r44_overnight_log.txt
    Write-Output "=== R44 $tag done, elapsed $((Get-Date) - $t0) ==="
}
Write-Output "=== R44 all cells done, total $((Get-Date) - $t0) ==="
