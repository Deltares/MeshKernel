@echo off

REM update_version.cmd replaces the depricated make_revision.bat
REM
REM All occurrences of make_revision.bat should be replaced by update_version.cmd
REM Argument conversion:
REM make_revision.bat m1 m2 m3 m4 m5
REM must be replaced by:
REM update_version.cmd m5 m2 m3
REM See issue Delft3D-16419



setlocal enableextensions

rem Program will replace %1 by the %1.svn and will replace VERSION_BUILD_NUMBER by a corresponding svn version using svnversion command
rem
rem %1 - path to the target source file
rem %2 - path to the folder to be used to check svnversion
rem %3 - Single file with version number information version_number.ini
rem %4 - --onlyifmissing: only regenerate when target source file does not exist (optional, default: off)

echo Generating version number in '%1' ...

set SCRIPT_DIRECTORY=%~dp0

set SVNVERSION="%SCRIPT_DIRECTORY%\..\..\third_party_open\subversion\bin\win32\svnversion.exe"
REM Temporariry fix until TeamCity is compatible with svn 1.8
set SVNVERSION17="%SCRIPT_DIRECTORY%\..\..\third_party_open\subversion\bin\win32-17\svnversion.exe"
set VN="%SCRIPT_DIRECTORY%..\..\third_party_open\version_number\bin\win32\version_number.exe"

echo version number %VN%
echo SVNVERSION %SVNVERSION%

set version=000000

rem Obtain the svn version number 
for /f "tokens=*" %%a in ('%SVNVERSION% %2') do set version=%%a
IF "%version%" == "000000" (
   REM use old svn version 1.7 for use on TeamCity Buildserver which does not support 1.8 yet
   for /f "tokens=*" %%a in ('%SVNVERSION17% %2') do set version=%%a
   )

rem In case version="123:128", use "128":
set version=%version:*:=%

rem ==========================================================================
rem If the source has been obtained using a svn export command, the "Unversioned directory"
rem string has been generated, but this cannot be used within *.rc files
rem Replace it using 000000 (only necessary on Windows systems)
rem ==========================================================================

IF "%version:~0,8%" == "exported" (
   set version=000000
)
IF "%version:~0,11%" == "Unversioned" (
   set version=000000
)
echo %0: %version%

IF "%4" == "--onlyifmissing" (
   IF EXIST "%1" (
      echo %0: Leaving existing file '%1' as is.
      goto end
   ) ELSE (
      echo %0: Create missing file '%1'.
   )
) ELSE (
   echo %0: Regenerating existing file '%1'.
)


if exist %1 (
	del %1
)

echo Arrived here %3
echo Arrived here "%1.svn"
echo Arrived here "%1"
rem Generate version number source module using version_number.exe
%VN% "%version%" "%3" "%1.svn" "%1"

:end
