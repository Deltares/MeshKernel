$ErrorActionPreference = 'Continue'

$results = [System.Collections.Generic.List[object]]::new()

# Run native programs through System.Diagnostics.Process. Windows PowerShell 5.1
# converts native stderr into ErrorRecords; capturing those with 2>&1 adds
# misleading NativeCommandError / source-line diagnostics to otherwise valid output.
function Invoke-NativeCommand {
    param(
        [Parameter(Mandatory = $true)][string]$FilePath,
        [Parameter(Mandatory = $true)][string[]]$Arguments
    )

    $quotedArguments = foreach ($argument in $Arguments) {
        '"' + ($argument -replace '"', '\\"') + '"'
    }

    $startInfo = New-Object System.Diagnostics.ProcessStartInfo
    $startInfo.FileName = $FilePath
    $startInfo.Arguments = $quotedArguments -join ' '
    $startInfo.UseShellExecute = $false
    $startInfo.RedirectStandardOutput = $true
    $startInfo.RedirectStandardError = $true
    $startInfo.CreateNoWindow = $true

    $process = New-Object System.Diagnostics.Process
    $process.StartInfo = $startInfo

    if (!$process.Start()) {
        throw "Could not start native command: $FilePath"
    }

    # Read both streams asynchronously to avoid blocking if either buffer fills.
    $stdoutTask = $process.StandardOutput.ReadToEndAsync()
    $stderrTask = $process.StandardError.ReadToEndAsync()
    $process.WaitForExit()

    [PSCustomObject]@{
        ExitCode = $process.ExitCode
        StdOut   = $stdoutTask.Result
        StdErr   = $stderrTask.Result
    }
}

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

function Test-SkippedCheck {
    param(
        [string]$Name,
        [string]$Reason
    )

    Write-Host "$Name : SKIPPED - $Reason" -ForegroundColor Yellow

    $results.Add([PSCustomObject]@{
        Check  = $Name
        Status = 'SKIP'
        Detail = $Reason
    })
}

# Find the repository root whether this script is launched from the
# repository root or from its conan subdirectory.
$repoRoot = $null
$searchStarts = @($PSScriptRoot, (Get-Location).Path) | Where-Object { $_ }

foreach ($start in $searchStarts) {
    $directory = Get-Item -LiteralPath $start -ErrorAction SilentlyContinue

    while ($directory -and !$repoRoot) {
        $candidate = $directory.FullName

        if ((Test-Path (Join-Path $candidate 'conanfile.py')) -and
            (Test-Path (Join-Path $candidate 'conan\profiles\windows-msvc-2022'))) {
            $repoRoot = $candidate
            break
        }

        $directory = $directory.Parent
    }

    if ($repoRoot) {
        break
    }
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
        ProcessArchitecture = if ([Environment]::Is64BitProcess) { 'x64' } else { 'x86' }
        Is64BitProcess      = [Environment]::Is64BitProcess
        ComputerName        = $env:COMPUTERNAME
    } | Format-List | Out-String

    if (![Environment]::Is64BitProcess) {
        throw 'Expected a 64-bit PowerShell process.'
    }
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

    if ($env:VSCMD_ARG_HOST_ARCH -ne 'x64') {
        throw 'Expected x64 Visual Studio host architecture.'
    }

    if ($env:VSCMD_ARG_TGT_ARCH -ne 'x64') {
        throw 'Expected x64 Visual Studio target architecture.'
    }
}

Test-Check 'Visual Studio 2022 Build Tools installation' {
    $expectedVsRoot = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\2022\BuildTools'
    $actualVsRoot = if ($env:VSINSTALLDIR) { $env:VSINSTALLDIR.TrimEnd('\') } else { $null }

    Write-Host "Expected: $expectedVsRoot"
    Write-Host "Active:   $actualVsRoot"

    if (!(Test-Path $expectedVsRoot)) {
        throw "Visual Studio 2022 Build Tools not found: $expectedVsRoot"
    }

    if (!$actualVsRoot -or
        ![string]::Equals($actualVsRoot, $expectedVsRoot, [StringComparison]::OrdinalIgnoreCase)) {
        throw "The active Developer environment is not Visual Studio 2022 Build Tools."
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

    Write-Host "`n--- Compiler and linker versions ---"
    $compilerVersion = (Get-Item -LiteralPath $compiler.Source).VersionInfo
    $linkerVersion = (Get-Item -LiteralPath $linker.Source).VersionInfo

    Write-Host "Compiler file version: $($compilerVersion.FileVersion)"
    Write-Host "Compiler product version: $($compilerVersion.ProductVersion)"
    Write-Host "Linker file version: $($linkerVersion.FileVersion)"
    Write-Host "Linker product version: $($linkerVersion.ProductVersion)"

    # Read version metadata instead of invoking cl.exe /Bv. In Windows
    # PowerShell 5.1, /Bv writes diagnostic text to stderr and exits with
    # D8003 when no source file is supplied; that is not a useful test here.
    if ($compilerVersion.FileVersion -notmatch '^19\.44\.') {
        throw "Expected MSVC compiler 19.44.x; found $($compilerVersion.FileVersion)."
    }

    if ($linkerVersion.FileVersion -notmatch '^14\.44\.') {
        throw "Expected MSVC linker 14.44.x; found $($linkerVersion.FileVersion)."
    }

    Write-Host "`n--- Compile and link ---"
    & $compiler.Source /nologo /EHsc /std:c++20 $source "/Fe:$exe"

    if ($LASTEXITCODE -ne 0) {
        throw "Compilation/linking failed with exit code $LASTEXITCODE"
    }

    if (!(Test-Path $exe)) {
        throw "Executable was not created: $exe"
    }

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

    if ($env:WindowsSDKVersion.TrimEnd('\') -ne '10.0.26100.0') {
        throw "Expected Windows SDK 10.0.26100.0, found $env:WindowsSDKVersion"
    }

    $expectedKernel32 = Join-Path $sdkRoot 'Lib\10.0.26100.0\um\x64\kernel32.lib'
    if (!(Test-Path $expectedKernel32)) {
        throw "Windows SDK x64 kernel32.lib not found: $expectedKernel32"
    }

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

    $rcVersion = (Get-Item -LiteralPath $expected).VersionInfo
    Write-Host "Resource compiler file version: $($rcVersion.FileVersion)"
}

Test-Check 'Resource compiler PATH ordering' {
    $env:PATH -split ';' |
        Where-Object { $_ -like '*Windows Kits\10\bin*' } |
        ForEach-Object { Write-Host $_ }

    $sdkVersion = $env:WindowsSDKVersion.TrimEnd('\')
    $expectedDirectory = Join-Path $env:WindowsSdkDir "bin\$sdkVersion\x86"
    $firstSdkBin = $env:PATH -split ';' |
        Where-Object { $_ -like '*Windows Kits\10\bin*' } |
        Select-Object -First 1

    $first = (Get-Command rc.exe -ErrorAction Stop).Source
    Write-Host "First SDK bin directory: $firstSdkBin"
    Write-Host "Resolved rc.exe: $first"

    if (!$firstSdkBin -or
        ![string]::Equals($firstSdkBin.TrimEnd('\'), $expectedDirectory.TrimEnd('\'), [StringComparison]::OrdinalIgnoreCase)) {
        throw "Expected the x86 SDK bin directory first in PATH: $expectedDirectory"
    }
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

Test-Check 'Conan cached default profile (informational)' {
    Write-Host 'Informational only: the MeshKernel build uses the explicit windows-msvc-2022 profile below.'
    conan profile show -pr default

    if ($LASTEXITCODE -ne 0) {
        throw "Conan returned exit code $LASTEXITCODE"
    }
}

Test-Check 'Conan Windows host and build profiles' {
    if (!$repoRoot) {
        throw 'MeshKernel repository root could not be located.'
    }

    $profile = Join-Path $repoRoot 'conan\profiles\windows-msvc-2022'
    $conan = Get-Command conan.exe -ErrorAction Stop
    $profileResult = Invoke-NativeCommand -FilePath $conan.Source -Arguments @(
        'profile', 'show', '-pr:h', $profile, '-pr:b', $profile
    )

    if ($profileResult.StdOut) {
        Write-Host $profileResult.StdOut.TrimEnd()
    }
    if ($profileResult.StdErr) {
        # Conan may write informational profile headings to stderr. Display
        # them as ordinary text rather than PowerShell NativeCommandError records.
        Write-Host $profileResult.StdErr.TrimEnd()
    }

    $profileOutput = $profileResult.StdOut + "`n" + $profileResult.StdErr

    if ($profileResult.ExitCode -ne 0) {
        throw "Conan profile show returned exit code $($profileResult.ExitCode)"
    }

    $expectedValues = @(
        'arch=x86_64',
        'build_type=Release',
        'compiler=msvc',
        'compiler.cppstd=20',
        'compiler.runtime=dynamic',
        'compiler.runtime_type=Release',
        'compiler.version=194',
        'os=Windows'
    )

    foreach ($value in $expectedValues) {
        if ($profileOutput -notmatch [regex]::Escape($value)) {
            throw "Expanded host/build profiles do not contain expected setting: $value"
        }
    }

    if ($profileOutput -notmatch 'cmake-msvc-rc-workaround\.cmake') {
        throw 'Expanded host/build profiles do not include the resource compiler workaround.'
    }

    Write-Host 'Host and build profiles both resolve with the expected Windows settings.'
}

# ------------------------------------------------------------------
# 7. MeshKernel and Conan project configuration
# ------------------------------------------------------------------

Test-Check 'MeshKernel repository files' {
    if (!$repoRoot) {
        throw 'Could not locate the MeshKernel repository root. Run this script from the repository or its conan directory.'
    }

    $requiredFiles = @(
        'conanfile.py',
        'conan\profiles\windows-msvc-2022',
        'conan\cmake\cmake-msvc-rc-workaround.cmake',
        'CMakeLists.txt'
    )

    Write-Host "Repository root: $repoRoot"

    $missing = @(
        foreach ($relativePath in $requiredFiles) {
            $fullPath = Join-Path $repoRoot $relativePath
            $exists = Test-Path $fullPath
            Write-Host "$exists : $fullPath"

            if (!$exists) {
                $fullPath
            }
        }
    )

    if ($missing.Count -gt 0) {
        throw "Missing repository files: $($missing -join '; ')"
    }
}

Test-Check 'Conan Windows profile settings and workaround' {
    if (!$repoRoot) {
        throw 'MeshKernel repository root could not be located.'
    }

    $profile = Join-Path $repoRoot 'conan\profiles\windows-msvc-2022'
    $profileText = Get-Content -LiteralPath $profile -Raw

    $requiredSettings = @(
        'arch=x86_64',
        'build_type=Release',
        'compiler=msvc',
        'compiler.cppstd=20',
        'compiler.runtime=dynamic',
        'compiler.runtime_type=Release',
        'compiler.version=194',
        'os=Windows'
    )

    foreach ($setting in $requiredSettings) {
        if ($profileText -notmatch ('(?m)^\s*' + [regex]::Escape($setting) + '\s*$')) {
            throw "Required Conan profile setting is missing: $setting"
        }
        Write-Host "Found: $setting"
    }

    if ($profileText -notmatch 'tools\.cmake\.cmaketoolchain:generator=Visual Studio 17 2022') {
        throw 'Conan profile does not select the Visual Studio 17 2022 generator.'
    }

    if ($profileText -notmatch 'tools\.cmake\.cmaketoolchain:user_toolchain=.*cmake-msvc-rc-workaround\.cmake') {
        throw 'Conan profile does not reference cmake-msvc-rc-workaround.cmake.'
    }

    $workaround = Join-Path $repoRoot 'conan\cmake\cmake-msvc-rc-workaround.cmake'
    $workaroundText = Get-Content -LiteralPath $workaround -Raw

    if ($workaroundText -notmatch 'CMAKE_VS_GLOBALS\s+"UseEnv=true"') {
        throw 'RC workaround does not set CMAKE_VS_GLOBALS to UseEnv=true.'
    }

    if ($workaroundText -notmatch 'Windows Kits/10/bin/10\.0\.26100\.0/x86/rc\.exe') {
        throw 'RC workaround does not select the x86 Windows SDK 10.0.26100.0 resource compiler.'
    }

    if ($workaroundText -notmatch 'CMAKE_TRY_COMPILE_PLATFORM_VARIABLES') {
        throw 'RC workaround does not propagate settings to try_compile() projects.'
    }

    Write-Host 'Visual Studio generator, profile, and RC workaround are configured.'
}

Test-Check 'CMake Conan provider can be disabled' {
    if (!$repoRoot) {
        throw 'MeshKernel repository root could not be located.'
    }

    $cmakeLists = Join-Path $repoRoot 'CMakeLists.txt'
    $cmakeText = Get-Content -LiteralPath $cmakeLists -Raw

    if ($cmakeText -notmatch 'option\s*\(\s*ENABLE_CONAN_PROVIDER') {
        throw 'CMakeLists.txt does not declare the ENABLE_CONAN_PROVIDER option.'
    }

    if ($cmakeText -notmatch 'if\s*\(\s*ENABLE_CONAN_PROVIDER\s*\)') {
        throw 'CMakeLists.txt does not conditionally enable the Conan provider.'
    }

    Write-Host 'ENABLE_CONAN_PROVIDER option and conditional include found.'
}

# The generated toolchain only exists after running `conan install`.
# Treat its absence as a skipped check so this script can also be used
# before the dependency-install step.
$conanToolchain = if ($repoRoot) {
    Join-Path $repoRoot 'build-conan-release\conan\conan_toolchain.cmake'
} else {
    $null
}

if ($conanToolchain -and (Test-Path $conanToolchain)) {
    Test-PathCheck 'Generated Conan CMake toolchain' $conanToolchain
}
else {
    Test-SkippedCheck 'Generated Conan CMake toolchain' 'Not generated yet; run the Conan install command first.'
}

# ------------------------------------------------------------------
# 8. Environment summary
# ------------------------------------------------------------------

Write-Host "`n========================================" -ForegroundColor Green
Write-Host ' Verification summary' -ForegroundColor Green
Write-Host '========================================' -ForegroundColor Green

$results | Format-Table -AutoSize

$passed = @($results | Where-Object Status -eq 'PASS')
$failed = @($results | Where-Object Status -eq 'FAIL')
$skipped = @($results | Where-Object Status -eq 'SKIP')

Write-Host "`nPassed:   $($passed.Count)"
Write-Host "Failed:   $($failed.Count)"
Write-Host "Skipped:  $($skipped.Count)"

if ($failed.Count -gt 0) {
    Write-Host "`nChecks requiring attention:" -ForegroundColor Red
    $failed | Format-Table -AutoSize
}
elseif ($skipped.Count -gt 0) {
    Write-Host "`nAll executed checks passed; some checks were skipped." -ForegroundColor Yellow
}
else {
    Write-Host "`nAll checks passed." -ForegroundColor Green
}
