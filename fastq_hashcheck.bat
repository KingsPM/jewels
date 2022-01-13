@echo off

REM Generates a hashfile for all files of a specified file extention saves to the same directory.
REM Written by Phil Davidson to enable hashchecking on fastq files by running on sequencer

set hashtype=md5
set file_type=fastq.gz
set output_file=fastq.md5

echo.
echo This script performs a hashcheck of %file_type% files and writes them to a file called %output_file%.
echo This %output_file% should then be copied to Snappy along with the %file_type% files.
echo.
echo Please provide the full fastq folder path i.e.
echo     On sequencers with Win10: d:\output\xxx_xxx_xxx\Alignment_1\xxx_xxx\Fastq
echo     On sequencers with Win7:  d:\Illumina\MiSeqOutput\xxx_xxx_xxx_xxx-xxx\Data\Intensities\BaseCalls
echo.
set /p source="Please enter fastq folder path: "

setlocal ENABLEDELAYEDEXPANSION

CALL :CreateHashfile !source!
pause
EXIT /B %ERRORLEVEL%


:: Functions
:CreateHashfile
echo.
echo Removing existing fastq.md5 if found...
del /q %~1\!output_file!
echo.

set total=0
for %%A in (%~1\*.!file_type!) do set /a total+=1
echo %total% !file_type!/s found in folder !source!
echo.

echo # Hashfile created on %date%%time% from hostname %computername% on !file_type! found in %~1 >> %~1\!output_file!
echo # Batchfile source code: https://git.kingspm.uk/KingsPM/jewels/src/fastq_hashcheck.bat >> %~1\!output_file!
echo # Author: "KCHBioinformatics <kch-tr.KCHBioinformatics@nhs.net>"  >> %~1\!output_file!

set /a count = 0
for /f "tokens=*" %%a in ('dir %~1\*.!file_type! /a:-d-l-h /b /s') do (
    set filepath=%%a
	for /r %%a in (!filepath!) do @set filename=%%~nxa
    for /f "eol=C tokens=* skip=1" %%b in ('certutil -hashfile "!filepath!" !hashtype!') do (
        set /a count += 1
        echo %%b !filename! ^(!count!/%total%^)
        echo %%b !filename! >> %~1\!output_file!
    )
)
echo.
echo !count! written to !output_file! file
echo.
EXIT /B 0
