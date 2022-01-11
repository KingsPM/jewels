@echo off

REM Generates a hashfile for all files of a specified file extention saves to the same directory.
REM Written by Phil Davidson to enable hashchecking on fastq files by running on sequencer
 
set hashtype=md5
set file_type=fastq.gz
set output_file=fastq.md5
set /p source="Directory Path: "

setlocal ENABLEDELAYEDEXPANSION

CALL :CreateHashfile !source!
pause
EXIT /B %ERRORLEVEL%


:: Functions
:CreateHashfile 
echo.
del /q %~1\!output_file!

echo # Hashfile created on %date%%time% from hostname %computername% >> %~1\!output_file!
echo # Batchfile source code: https://git.kingspm.uk/KingsPM/jewels/src/branch/develop/fastq_md5_hashcheck.bat >> %~1\!output_file!
echo # Author: "KCHBioinformatics <kch-tr.KCHBioinformatics@nhs.net>"  >> %~1\!output_file!

set /a count = 0
for /f "tokens=*" %%a in ('dir %~1\*.!file_type! /a:-d-l-h /b /s') do (
    set filepath=%%a
	for /r %%a in (!filepath!) do @set filename=%%~nxa
    for /f "eol=C tokens=* skip=1" %%b in ('certutil -hashfile "!filepath!" !hashtype!') do (
        set /a count += 1
        echo %%b !filename!
        echo %%b !filename! >> %~1\!output_file!
    )
)
echo.
echo !count! files found and written to !output_file! file
echo.
EXIT /B 0
