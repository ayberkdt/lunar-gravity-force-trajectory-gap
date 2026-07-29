param(
    [Parameter(Mandatory=$true)]
    [int]$QueuePid
)

$pythonExe = 'D:\Masaustu\LUNAR_SIMULATION\.venv\Scripts\python.exe'
$script = 'C:\Users\ayber\Desktop\Makale\codebase\python_codes\rev10_queue_postprocess.py'
Write-Output "[postprocess-watcher] waiting for queue PID $QueuePid"
Wait-Process -Id $QueuePid -ErrorAction SilentlyContinue
Write-Output '[postprocess-watcher] queue ended; generating aggregate'
& $pythonExe -u $script
Write-Output "[postprocess-watcher] exit=$LASTEXITCODE"
