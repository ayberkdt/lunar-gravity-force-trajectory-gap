@echo off
rem R50 campaign launcher, run through Task Scheduler.
rem
rem Two launches died before this file existed: Start-Process from the agent
rem shell was killed with the session (DBG_TERMINATE_PROCESS) and a
rem Win32_Process.Create launch took a console control event
rem (STATUS_CONTROL_C_EXIT). A scheduled task is owned by the scheduler
rem service, so no console event and no session teardown reaches it.
cd /d "C:\Users\ayber\Desktop\Makale\codebase\python_codes"
"C:\Users\ayber\AppData\Local\Programs\Python\Python312\python.exe" rev50_campaign.py --stop-at 2026-08-13T02:00:00 --workers 11 >> "C:\Users\ayber\Desktop\Makale\codebase\output\r50_campaign_20260812c.stdout.log" 2>&1
