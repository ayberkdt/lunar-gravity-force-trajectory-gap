# Overnight chain v2, relaunched ~23:30 local 2026-08-09. Cutoffs local.
# Stage 1: R44 equal-work re-match, tighter level, 8 cells (resume-safe).
# Stage 2: O48 interior measured-time panel, SERIAL (machine idle inside the
#          chain), cutoff 06:45.
# Stage 3: R14 design-B beta=2.00 then 3.00 if time remains (skip past 07:00;
#          rev14 --deadline is ISO and naive = UTC: 05:15 UTC = 08:15 local).
$ErrorActionPreference = "Continue"
Set-Location $PSScriptRoot
$log = "..\metrics\r44_overnight_log.txt"
$t0 = Get-Date
"=== OVERNIGHT v2 start $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') ===" |
    Tee-Object -Append -FilePath $log

# ---- Stage 1: R44 ----
$cells = @(
    @{d = "A"; b = "1.00"}, @{d = "B"; b = "1.00"},
    @{d = "A"; b = "0.75"}, @{d = "B"; b = "0.75"},
    @{d = "A"; b = "0.50"}, @{d = "B"; b = "0.50"},
    @{d = "A"; b = "1.25"}, @{d = "B"; b = "1.50"}
)
foreach ($c in $cells) {
    "=== R44 $($c.d)@$($c.b) start $(Get-Date -Format HH:mm:ss) ===" |
        Tee-Object -Append -FilePath $log
    python rev44_equal_work_tighter.py run --design $c.d --beta $c.b `
        --workers 10 --deadline-min 60 2>&1 |
        Tee-Object -Append -FilePath $log
}
"=== R44 done, elapsed $((Get-Date) - $t0) ===" |
    Tee-Object -Append -FilePath $log

# ---- Stage 2: O48 serial measured-time panel ----
"=== O48 pipeline start $(Get-Date -Format HH:mm:ss) ===" |
    Tee-Object -Append -FilePath $log
python rev48_interior_timing.py pipeline --cutoff 2026-08-10T06:45 2>&1 |
    Tee-Object -Append -FilePath $log
"=== O48 done, elapsed $((Get-Date) - $t0) ===" |
    Tee-Object -Append -FilePath $log

# ---- Stage 3: R14 design-B upper-budget cells ----
foreach ($b in @("2.00", "3.00")) {
    if ((Get-Date) -gt (Get-Date "2026-08-10 07:00")) {
        "=== SKIP R14 B@$b (past 07:00 local) ===" |
            Tee-Object -Append -FilePath $log
        continue
    }
    "=== R14 B@$b start $(Get-Date -Format HH:mm:ss) ===" |
        Tee-Object -Append -FilePath $log
    python rev14_budget_trajectory.py run --design B --beta $b `
        --workers 10 --deadline 2026-08-10T05:15 2>&1 |
        Tee-Object -Append -FilePath $log
}
"=== OVERNIGHT v2 done, total $((Get-Date) - $t0) ===" |
    Tee-Object -Append -FilePath $log
