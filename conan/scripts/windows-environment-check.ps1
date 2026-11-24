
$ErrorActionPreference = 'Continue'

$results = [System.Collections.Generic.List[object]]::new()

function Test-Check {
    param(
        [string]$Name,
        [scriptblock]$Action
    )

    Write-Host "`n=== $Name ===" -ForegroundColor Cyan

    try {
        # Reset the native exit code before each independent check.
        $global:LASTEXITCODE = 0

        $output = & $Action 2>&1
        $exitCode = $global:LASTEXITCODE

        foreach ($line in $output) {
            Write-Host $line
        }

        $ok = ($exitCode -eq 0)

        $results.Add([PSCustomObject]@{
            Check  = $Name
            Status = if ($ok) { 'PASS' } else { 'FAIL' }
            Detail = "Exit code: $exitCode"
        })
    }
    catch {
        Write-Host $_ -ForegroundColor Red

        $results.Add([PSCustomObject]@{
            Check  = $Name
            Status = 'FAIL'
            Detail = $_.Exception.Message
        })
    }
}

function Test-PathCheck {
    param(
        [string]$Name,
        [string]$Path
    )

    $exists = Test-Path $Path
    Write-Host "$Name : $Path"
    Write-Host "Exists: $exists"

    $results.Add([PSCustomObject]@{
        Check  = $Name
        Status = if ($exists) { 'PASS' } else { 'FAIL' }
        Detail = $Path
    })
}

Write-Host '========================================' -ForegroundColor Green
Write-Host ' MeshKernel Windows build environment' -ForegroundColor Green
Write-Host '========================================' -ForegroundColor Green

# ------------------------------------------------------------------
# 1. Operating system and process architecture
# ------------------------------------------------------------------

Test-Check 'Operating system' {
    Get-CimInstance Win32_OperatingSystem |
        Select-Object Caption, Version, BuildNumber |
        Format-List | Out-String
}

Test-Check 'PowerShell process architecture' {
    [PSCustomObject]@{
        PowerShellVersion   = $PSVersionTable.PSVersion.ToString()
        ProcessArchitecture = [Runtime.InteropServices.RuntimeInformation]::ProcessArchitecture
        Is64BitProcess      = [Environment]::Is64BitProcess
        ComputerName        = $env:COMPUTERNAME
    } | Format-List | Out-String
}

# ------------------------------------------------------------------
# 2. Visual Studio Developer environment
# ------------------------------------------------------------------

Test-Check 'Visual Studio environment variables' {
    Get-ChildItem Env: |
        Where-Object {
            $_.Name -match '^(VS|VC|VSCMD|WindowsSdk|UniversalCRTSdk|UCRTVersion|INCLUDE$|LIB$|LIBPATH$)'
        } |
        Sort-Object Name |
        Format-Table -AutoSize | Out-String
}

Test-Check 'Visual Studio installation' {
    $vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"

    if (!(Test-Path $vswhere)) {
        throw "vswhere.exe not found: $vswhere"
    }

    & $vswhere -latest -products '*' -property installationPath

    if ($LASTEXITCODE -ne 0) {
        throw "vswhere.exe returned exit code $LASTEXITCODE"
    }
}

Test-Check 'Developer command prompt configuration' {
    [PSCustomObject]@{
        VSCMD_VER           = $env:VSCMD_VER
        VSCMD_ARG_HOST_ARCH = $env:VSCMD_ARG_HOST_ARCH
        VSCMD_ARG_TGT_ARCH  = $env:VSCMD_ARG_TGT_ARCH
        VSINSTALLDIR        = $env:VSINSTALLDIR
        VCINSTALLDIR        = $env:VCINSTALLDIR
        WindowsSdkDir       = $env:WindowsSdkDir
        WindowsSDKVersion   = $env:WindowsSDKVersion
    } | Format-List | Out-String

    if (!$env:VSCMD_VER -or !$env:VSINSTALLDIR -or !$env:VCINSTALLDIR) {
        throw 'Visual Studio Developer environment is not fully initialised.'
    }

    if ($env:VSCMD_ARG_TGT_ARCH -ne 'x64') {
        throw 'Expected x64 Visual Studio target architecture.'
    }
}

# ------------------------------------------------------------------
# 3. MSVC compiler and linker
# ------------------------------------------------------------------

Test-Check 'MSVC compiler location' {
    $compiler = Get-Command cl.exe -ErrorAction Stop
    Write-Host $compiler.Source
}

Test-Check 'MSVC linker location' {
    $linker = Get-Command link.exe -ErrorAction Stop
    Write-Host $linker.Source
}

Test-Check 'MSVC compiler and linker version' {
    $compiler = Get-Command cl.exe -ErrorAction Stop
    $linker = Get-Command link.exe -ErrorAction Stop

    Write-Host "Compiler: $($compiler.Source)"
    Write-Host "Linker:   $($linker.Source)"

    # /Bv prints the compiler version but returns exit code 2
    # without a source file. Use an actual compile/link test instead.

    $temp = Join-Path $env:TEMP 'meshkernel-msvc-check'
    New-Item -ItemType Directory -Force -Path $temp | Out-Null

    $source = Join-Path $temp 'test.cpp'
    $exe = Join-Path $temp 'test.exe'

    @'
#include <Windows.h>

int main()
{
    return GetCurrentProcessId() == 0;
}
'@ | Set-Content -Path $source -Encoding ASCII

    Write-Host "`n--- Compiler version ---"
    & $compiler.Source /nologo /Bv

    # Ignore the expected nonzero exit from /Bv without a source file.
    $global:LASTEXITCODE = 0

    Write-Host "`n--- Compile and link ---"
    & $compiler.Source /nologo /EHsc /std:c++20 $source "/Fe:$exe"

    if ($LASTEXITCODE -ne 0) {
        throw "Compilation/linking failed with exit code $LASTEXITCODE"
    }

    if (!(Test-Path $exe)) {
        throw "Executable was not created: $exe"
    }

    Write-Host "`n--- Linker version ---"
    $linkerVersion = (Get-Item $linker.Source).VersionInfo
    Write-Host "File version: $($linkerVersion.FileVersion)"
    Write-Host "Product version: $($linkerVersion.ProductVersion)"

    Write-Host "`n--- Execute test program ---"
    & $exe

    if ($LASTEXITCODE -ne 0) {
        throw "Test executable returned exit code $LASTEXITCODE"
    }

    Write-Host 'Compile, link and execution successful.'
}

Test-Check 'MSBuild location and version' {
    $msbuild = Get-Command MSBuild.exe -ErrorAction Stop

    Write-Host $msbuild.Source
    & $msbuild.Source -version -nologo

    if ($LASTEXITCODE -ne 0) {
        throw "MSBuild returned exit code $LASTEXITCODE"
    }
}

# ------------------------------------------------------------------
# 4. Windows SDK and resource compiler
# ------------------------------------------------------------------

Test-Check 'Windows SDK installation' {
    $sdkRoot = $env:WindowsSdkDir.TrimEnd('\')

    if (!(Test-Path $sdkRoot)) {
        throw "Windows SDK root not found: $sdkRoot"
    }

    Write-Host "SDK root: $sdkRoot"
    Write-Host "SDK version: $env:WindowsSDKVersion"

    Get-ChildItem "$sdkRoot\bin" -Directory |
        Select-Object Name, FullName |
        Format-Table -AutoSize | Out-String
}

Test-Check 'Resource compiler selection' {
    $rc = Get-Command rc.exe -ErrorAction Stop
    Write-Host "Selected rc.exe: $($rc.Source)"

    $sdkVersion = $env:WindowsSDKVersion.TrimEnd('\')
    $expected = Join-Path $env:WindowsSdkDir "bin\$sdkVersion\x86\rc.exe"

    Write-Host "Expected rc.exe: $expected"

    if ($rc.Source -ne $expected) {
        throw "Wrong resource compiler selected. Expected $expected"
    }

    if (!(Test-Path $expected)) {
        throw "Resource compiler does not exist: $expected"
    }

    & $expected /?

    # rc.exe /? is a usage query; do not use its exit code
    # as an installation or functionality test.
}

Test-Check 'Resource compiler PATH ordering' {
    $env:PATH -split ';' |
        Where-Object { $_ -like '*Windows Kits\10\bin*' } |
        ForEach-Object { Write-Host $_ }

    $first = (Get-Command rc.exe -ErrorAction Stop).Source
    Write-Host "Resolved rc.exe: $first"
}

Test-Check 'Windows SDK libraries and headers' {
    $sdkVersion = $env:WindowsSDKVersion.TrimEnd('\')
    $sdkRoot = $env:WindowsSdkDir.TrimEnd('\')

    $includeRoot = Join-Path $sdkRoot "Include\$sdkVersion"
    $libRoot = Join-Path $sdkRoot "Lib\$sdkVersion"
    $binRoot = Join-Path $sdkRoot "bin\$sdkVersion"

    $paths = @(
        "$includeRoot\um\Windows.h",
        "$includeRoot\shared\winapifamily.h",
        "$libRoot\um\x64\kernel32.lib",
        "$libRoot\ucrt\x64\ucrt.lib",
        "$binRoot\x86\rc.exe"
    )

    $missing = @(
        foreach ($path in $paths) {
            $exists = Test-Path $path
            Write-Host "$exists : $path"

            if (!$exists) {
                $path
            }
        }
    )

    if ($missing.Count -gt 0) {
        throw "Missing Windows SDK files: $($missing -join '; ')"
    }
}

# ------------------------------------------------------------------
# 5. Include and library paths
# ------------------------------------------------------------------

Test-Check 'Compiler INCLUDE paths' {
    $env:INCLUDE -split ';' |
        ForEach-Object { Write-Host $_ }

    if (!$env:INCLUDE) {
        throw 'INCLUDE is not set.'
    }
}

Test-Check 'Compiler LIB paths' {
    $env:LIB -split ';' |
        ForEach-Object { Write-Host $_ }

    if (!$env:LIB) {
        throw 'LIB is not set.'
    }
}

# ------------------------------------------------------------------
# 6. Build tools
# ------------------------------------------------------------------

Test-Check 'CMake' {
    cmake --version

    if ($LASTEXITCODE -ne 0) {
        throw "CMake returned exit code $LASTEXITCODE"
    }
}

Test-Check 'Git' {
    git --version

    if ($LASTEXITCODE -ne 0) {
        throw "Git returned exit code $LASTEXITCODE"
    }
}

Test-Check 'Python' {
    python --version

    if ($LASTEXITCODE -ne 0) {
        throw "Python returned exit code $LASTEXITCODE"
    }

    python -c "import sys; print(sys.executable)"

    if ($LASTEXITCODE -ne 0) {
        throw "Python returned exit code $LASTEXITCODE"
    }
}

Test-Check 'Conan' {
    conan --version

    if ($LASTEXITCODE -ne 0) {
        throw "Conan returned exit code $LASTEXITCODE"
    }
}

Test-Check 'Conan profiles' {
    conan profile list

    if ($LASTEXITCODE -ne 0) {
        throw "Conan returned exit code $LASTEXITCODE"
    }
}

Test-Check 'Conan default profile' {
    conan profile show -pr default

    if ($LASTEXITCODE -ne 0) {
        throw "Conan returned exit code $LASTEXITCODE"
    }
}

# ------------------------------------------------------------------
# 7. Environment summary
# ------------------------------------------------------------------

Write-Host "`n========================================" -ForegroundColor Green
Write-Host ' Verification summary' -ForegroundColor Green
Write-Host '========================================' -ForegroundColor Green

$results | Format-Table -AutoSize

$failed = @($results | Where-Object Status -eq 'FAIL')

Write-Host "`nPassed: $($results.Count - $failed.Count)"
Write-Host "Failed: $($failed.Count)"

if ($failed.Count -gt 0) {
    Write-Host "`nChecks requiring attention:" -ForegroundColor Red
    $failed | Format-Table -AutoSize
}
else {
    Write-Host "`nAll checks passed." -ForegroundColor Green
}
