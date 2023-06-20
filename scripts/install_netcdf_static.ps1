
<#
.Description
This script builds and installs the HDF5 and NetCDF static libraries.
These libraries are required for building the MeshKernel unit tests.
.PARAMETER BuildType
Sets the build type.
.PARAMETER WorkDir
Path to the directory where packages are downloaded and extracted, repositories are cloned, and targets are built.
.PARAMETER InstallDir
Path to the directory where the libraries are installed.
.PARAMETER GitTags
Git tags to checkout the desired git branches.
Care must be taken when setting the tag strings. Different repositories use different tag formats.
Let X, Y and Z be the major, minor and patch release numbers, respectively. The tags are formatted as follows:
zlib     : 'vX.Y.Z'
curl     : 'curl-X_Y_Z'
hdf5     : 'hdf5-X_Y_Z'
netcdf_c : 'vX.Y.Z'
For valid tags check:
zlib     : https://github.com/madler/zlib/tags
curl     : https://github.com/curl/curl/tags
hdf5     : https://github.com/HDFGroup/hdf5/tags
netcdf_c : https://github.com/Unidata/netcdf-c/tags
.PARAMETER ParallelJobs
The maximum number of concurrent processes to use when building.
.PARAMETER Clean
Removes the work directory upon finishing the installation.
.SYNOPSIS
Used to build the HDF5 and NetCDF static libraries.
#>

Param(
    [Parameter(Mandatory = $true)] [string] $WorkDir,
    [Parameter(Mandatory = $true)] [string] $InstallDir,
    [Parameter(Mandatory = $false)] [ValidateSet('Release', 'Debug', 'RelWithDebInfo')] [string] $BuildType = 'Release',
    [Parameter(Mandatory = $false)] [hashtable]$GitTags = @{ `
            zlib     = 'v1.2.13'; `
            curl     = 'curl-7_88_1'; `
            hdf5     = 'hdf5-1_14_0'; `
            netcdf_c = 'v4.9.1'
    },
    [Parameter(Mandatory = $false)] [ValidateRange(1, [int]::MaxValue)] [int] $ParallelJobs = 6,
    [Parameter(Mandatory = $false)] [Switch] $Clean = $False
)

#$ErrorActionPreference = "Stop"

if (-not(Test-Path -Path $WorkDir)) {
    New-Item $WorkDir -Type Directory
}
$WorkDir = Resolve-Path $WorkDir

$TransriptPath = (Join-Path $WorkDir 'log.txt')
Start-Transcript -Path $TransriptPath

Function Invoke-Terminate {
    Stop-Transcript
}

# Check if OS is supported
if ( -not ($IsLinux -or $IsMacOS -or $IsWindows)) {
    Write-Error ('Unsupported operating system. The follwoing are supported: Linux, macOS and Windows.')
    Invoke-Terminate
    Exit
}

# Check version, powershell 7+ is supported
$MajorPSVers = $PSVersionTable.PSVersion.Major
if ([int]$MajorPSVers -lt 7) {
    Write-Error ('Found PowerShell version $MajorPSVers. Version 7+ is required.')
    Invoke-Terminate
    Exit
}

if (-not(Test-Path -Path $InstallDir)) {
    New-Item $InstallDir -Type Directory
}
$InstallDir = Resolve-Path $InstallDir

Write-Host 'Work directory              : ' $WorkDir
Write-Host 'Installation directory      : ' $InstallDir
Write-Host 'Build type                  : ' $BuildType
Write-Host 'Parallel jobs per build     : ' $ParallelJobs
Write-Host 'Tagged branches to checkout : ' ($GitTags | Out-String)

# No need to download and extract m4 when platform is WIN32
if ($IsWindows) {
    # ----------------------------------------------------------------------------------------
    # Download
    # ----------------------------------------------------------------------------------------

    $DownloadDir = (Join-Path $WorkDir 'download')
    New-Item -Force $DownloadDir -Type Directory

    $WebClient = New-Object System.Net.WebClient

    # M4 URL
    $M4BaseURL = 'http://downloads.sourceforge.net/gnuwin32'

    # Binary URL
    $M4BinFileName = 'm4-1.4.14-1-bin.zip'
    $M4BinURL = (@($M4BaseURL; $M4BinFileName) -Join '/')
    $M4BinDownloadPath = (Join-Path $DownloadDir $M4BinFileName)
    if (-not(Test-Path -Path $M4BinDownloadPath)) {
        $WebClient.DownloadFile($M4BinURL, $M4BinDownloadPath)
    }

    # Dependencies URL 
    $M4DepFileName = 'm4-1.4.14-1-dep.zip'
    $M4DepURL = (@($M4BaseURL; $M4DepFileName) -Join '/')
    $M4DepDownloadPath = (Join-Path $DownloadDir $M4DepFileName)
    if (-not(Test-Path -Path $M4DepDownloadPath)) {
        $WebClient.DownloadFile($M4DepURL, $M4DepDownloadPath)
    }

    # ----------------------------------------------------------------------------------------
    # Extract
    # ----------------------------------------------------------------------------------------

    $ExtractDir = (Join-Path $WorkDir 'extract')
    New-Item -Force $ExtractDir -Type Directory

    # M4
    $M4BinExtractPath = (Join-Path $ExtractDir 'm4-1.4.14-1-bin')
    Expand-Archive -Force  $M4BinDownloadPath -DestinationPath $M4BinExtractPath
    $env:Path += (';' + (Join-Path $M4BinExtractPath 'bin'))
    $M4DepExtractPath = (Join-Path $ExtractDir 'm4-1.4.14-1-dep')
    Expand-Archive -Force  $M4DepDownloadPath -DestinationPath $M4DepExtractPath
    $env:Path += (';' + (Join-Path $M4DepExtractPath 'bin'))
}

# ----------------------------------------------------------------------------------------
# Clone
# ----------------------------------------------------------------------------------------

$ReposDir = (Join-Path $WorkDir 'repositories')
New-Item -Force $ReposDir -Type Directory

Function Invoke-CloneRepoAndCheckoutTag {
    Param
    (
        [Parameter(Mandatory = $true)] [string] $Repo,
        [Parameter(Mandatory = $true)] [string] $Tag,
        [Parameter(Mandatory = $true)] [string] $Destination
    )

    if (-not(Test-Path -Path $Destination)) {
        #Remove-Item -Recurse -Force $Destination
        New-Item $Destination -Type Directory
    
        # clone the repo to the specified destination
        $Clone = ('git clone {0} {1}' -f $Repo, $Destination)
        Write-Host $Clone
        Invoke-Expression $Clone
        $ExitCode = $LASTEXITCODE
        if ( $ExitCode -gt 0 ) {
            Write-Error ('Failed to clone {0} to {1}. Exit code = {2}.' `
                    -f $Repo, $Destination, $ExitCode.ToString())
            Invoke-Terminate
            Exit
        }
        # enter the repo and fetch all tags
        $FetchTags = ('git -C {0} fetch --all --tags' -f $Destination)
        Write-Host $FetchTags
        Invoke-Expression $FetchTags
        $ExitCode = $LASTEXITCODE
        if ( $ExitCode -gt 0 ) {
            Write-Error ('Failed to fetch all tags from {0}. Exit code = {1}.' `
                    -f $Repo, $ExitCode.ToString())
            Invoke-Terminate
            Exit
        }

        # check out the brang with the specified tag
        $CheckOutTaggedBranch = ('git -C {0} checkout tags/{1} -b {1}' -f $Destination, $Tag)
        Write-Host $CheckOutTaggedBranch
        Invoke-Expression $CheckOutTaggedBranch
        $ExitCode = $LASTEXITCODE
        if ( $ExitCode -gt 0 ) {
            Write-Error ('Failed to checkout tags/{0} from {1}. Exit code = {2}.' `
                    -f $Tag, $Repo, $ExitCode.ToString())
            Invoke-Terminate
            Exit
        }
    }
}

# ZLIB
$ZLIB = 'zlib'
$ZLIBRepo = 'https://github.com/madler/zlib.git'
$ZLIBSrcDir = (Join-Path $ReposDir $ZLIB)
Invoke-CloneRepoAndCheckoutTag -Repo $ZLIBRepo -Tag $GitTags.zlib -Destination $ZLIBSrcDir

#Curl
if ($IsWindows) {
    $Curl = 'curl'
    $CurlRepo = 'https://github.com/curl/curl.git'
    $CurlSrcDir = (Join-Path $ReposDir $Curl)
    Invoke-CloneRepoAndCheckoutTag -Repo $CurlRepo -Tag $GitTags.curl -Destination $CurlSrcDir
}

# HDF5
$HDF5 = 'hdf5'
$HDF5Repo = 'https://github.com/HDFGroup/hdf5.git'
$HDF5SrcDir = (Join-Path $ReposDir $HDF5)
Invoke-CloneRepoAndCheckoutTag -Repo $HDF5Repo -Tag $GitTags.hdf5 -Destination $HDF5SrcDir

# NetCDF
$NetCDF = 'netcdf-c'
$NetCDFRepo = 'https://github.com/Unidata/netcdf-c.git'
$NetCDFSrcDir = (Join-Path $ReposDir $NetCDF)
Invoke-CloneRepoAndCheckoutTag -Repo $NetCDFRepo -Tag $GitTags.netcdf_c -Destination $NetCDFSrcDir

# ----------------------------------------------------------------------------------------
# build
# ----------------------------------------------------------------------------------------
$BuildDir = (Join-Path $WorkDir 'build')
New-Item -Force $BuildDir -Type Directory

$LocalInstallDir = (Join-Path $WorkDir 'install')

# adapted from https://www.powershellgallery.com/packages/WebKitDev/0.1.5/Content/Functions%5CInvoke-CMakeBuild.ps1
Function Invoke-BuildAndInstall {
    Param(
        [Parameter(Mandatory = $true)] [string] $SrcDir,
        [Parameter(Mandatory = $true)] [string] $BuildDir,
        [Parameter(Mandatory = $true)] [string] $InstallDir,
        [Parameter(Mandatory = $true)] [ValidateSet('Release', 'Debug', 'RelWithDebInfo')] [string] $BuildType,
        [Parameter(Mandatory = $false)] [ValidateRange(1, [int]::MaxValue)] [int] $ParallelJobs,
        [Parameter(Mandatory = $false)] [string[]] $Options = @()
    )
    
    if (-not(Test-Path -Path $BuildDir)) {
        New-Item $BuildDir -Type Directory
    }

    if (-not(Test-Path -Path $InstallDir)) {
        New-Item $InstallDir -Type Directory
    }
  
    # Configure
    $ConfigArgs = @('-S {0}' -f $SrcDir)
    $ConfigArgs += ('-B {0}' -f $BuildDir)
    $ConfigArgs += ('-DCMAKE_INSTALL_PREFIX={0}' -f $InstallDir)
    $ConfigArgs += ('-DCMAKE_BUILD_TYPE={0}' -f $BuildType)
    $ConfigArgs += $Options
    $Configure = ('cmake {0}' -f ($ConfigArgs -Join ' '))
    Write-Host $Configure
    Invoke-Expression $Configure
    $ExitCode = $LASTEXITCODE
    if ( $ExitCode -gt 0 ) {
        Write-Error ('Failed to configure the build {0}. Exit code = {1}.' `
                -f $SrcDir, $ExitCode.ToString())
        Invoke-Terminate
        Exit
    }

    # Build
    $BuildArgs += @('--build', $BuildDir, '--config', $BuildType, '--parallel', $ParallelJobs.ToString())
    $Build = ('cmake {0}' -f ($BuildArgs -Join ' '))
    Write-Host $Build
    Invoke-Expression $Build
    $ExitCode = $LASTEXITCODE
    if ( $ExitCode -gt 0 ) {
        Write-Error ('Failed to build {0} in {1}. Exit code = {2}.' `
                -f $SrcDir, $BuildDir, $ExitCode.ToString())
        Invoke-Terminate
        Exit
    }

    # Install
    $InstalldArgs += @('--install', $BuildDir)
    $Install = ('cmake {0}' -f ($InstalldArgs -Join ' '))
    Write-Host $Install
    Invoke-Expression $Install
    $ExitCode = $LASTEXITCODE
    if ( $ExitCode -gt 0 ) {
        Write-Error ('Failed to install {0} in {1}. Exit code = {2}.' `
                -f $BuildDir, $InstallDir, $ExitCode.ToString())
        Exit
    }
}

# ZLIB
$ZLIBBuildDir = (Join-Path $BuildDir $ZLIB)
$ZLIBInstallDir = (Join-Path $LocalInstallDir $ZLIB)
Invoke-BuildAndInstall `
    -SrcDir $ZLIBSrcDir `
    -BuildDir $ZLIBBuildDir `
    -InstallDir $ZLIBInstallDir `
    -ParallelJobs $ParallelJobs `
    -BuildType $BuildType
$env:Path += (';' + $ZLIBInstallDir)

# Curl
if ($IsWindows) {
    $CurlBuildDir = (Join-Path $BuildDir $Curl)
    $CurlInstallDir = (Join-Path $LocalInstallDir $Curl)
    Invoke-BuildAndInstall `
        -SrcDir $CurlSrcDir `
        -BuildDir $CurlBuildDir `
        -InstallDir $CurlInstallDir `
        -ParallelJobs $ParallelJobs `
        -BuildType $BuildType
    $env:Path += (';' + $CurlInstallDir)
}

# HDF5
$HDF5BuildDir = (Join-Path $BuildDir $HDF5)
$HDF5InstallDir = (Join-Path $LocalInstallDir $HDF5)
$ZlibIncludeDir = (Join-Path $ZLIBInstallDir 'include')
if ($IsLinux -or $IsMacOS) {
    $ZlibStaticLibrary = (Join-Path (Join-Path $ZLIBInstallDir 'lib') 'libz.a')
}
else {
    $ZlibStaticLibrary = (Join-Path (Join-Path $ZLIBInstallDir 'lib') 'zlib.lib')
}

$HDF5CMakeBuildOptions = @( `
        '-DBUILD_STATIC_LIBS:BOOL=ON', `
        '-DBUILD_SHARED_LIBS:BOOL=OFF', `
        '-DHDF5_BUILD_TOOLS:BOOL=OFF', `
        '-DHDF5_BUILD_EXAMPLES:BOOL=OFF', `
        '-DBUILD_TESTING:BOOL=OFF', `
        '-DZLIB_USE_EXTERNAL:BOOL=OFF', `
        '-DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON', `
    ('-DZLIB_ROOT={0}' -f $ZLIBInstallDir), `
    ('-DZLIB_INCLUDE_DIR:PATH={0}' -f $ZlibIncludeDir), `
    ('-DZLIB_LIBRARY:FILEPATH={0}' -f $ZlibStaticLibrary) `
)

Invoke-BuildAndInstall `
    -SrcDir $HDF5SrcDir `
    -BuildDir $HDF5BuildDir `
    -InstallDir $HDF5InstallDir `
    -BuildType $BuildType `
    -ParallelJobs $ParallelJobs `
    -Options $HDF5CMakeBuildOptions

$env:Path += ('; ' + $HDF5InstallDir)
$env:HDF5_ROOT = $HDF5InstallDir

# NetCDF
$NetCDFBuildDir = (Join-Path $BuildDir $NetCDF)
$NetCDFInstallDir = (Join-Path $InstallDir $NetCDF)
$NetCDFCMakeBuildOptions = @(`
        '-DBUILD_SHARED_LIBS=OFF', `
        '-DENABLE_NETCDF_4=ON', `
        '-DENABLE_DAP=OFF', `
        '-DENABLE_BYTERANGE=OFF', `
    ('-DZLIB_ROOT={0}' -f $ZLIBInstallDir) `
)
if ($IsLinux -or $IsMacOS) {
    $NetCDFCMakeBuildOptions += '-DCMAKE_POSITION_INDEPENDENT_CODE=ON'
}

Invoke-BuildAndInstall `
    -SrcDir $NetCDFSrcDir `
    -BuildDir $NetCDFBuildDir `
    -InstallDir $NetCDFInstallDir `
    -BuildType $BuildType `
    -ParallelJobs $ParallelJobs `
    -Options $NetCDFCMakeBuildOptions

# Some post-build manual installations... and arguably a terrible idea. This might break in future HDF5 or NetCDF releases.
# 
# Under WIN32, NetCDF links with -lhdf5-static -lhdf5_hl-static -lzlib. See:
# $NetCDFInstallDir/lib/libnetcdf.settings and $NetCDFInstallDir/lib/pkgconfig/netxdf.pc
# So we copy the static ZLIB and HDF5 lib dependnecies to $NetCDFInstallDir and rename them accordingly,
# and finally edit the list of public interface libararies in $NetCDFInstallDir/lib/cmake/netCDF/netCDFTargets.cmake.
# Under Linux, we simply copy static and shared libs without renaming.
Function Invoke-Post-Build-Steps() {
    # Copy all necessary static libraries from the local instalaltion directory to the netcdf lib dir
    $NetCDFLibDir = (Join-Path $NetCDFInstallDir 'lib')
    if ($IsLinux -or $IsMacOS) {
        $NetCDFBinDir = (Join-Path $NetCDFInstallDir 'bin')
        Copy-Item (Join-Path $ZLIBInstallDir 'lib' 'libz.a')       -Destination (Join-Path $NetCDFLibDir 'libz.a')
        Copy-Item (Join-Path $HDF5InstallDir 'lib' 'libhdf5.a')    -Destination (Join-Path $NetCDFLibDir 'libhdf5.a')
        Copy-Item (Join-Path $HDF5InstallDir 'lib' 'libhdf5_hl.a') -Destination (Join-Path $NetCDFLibDir 'libhdf5_hl.a')
        if ($IsLinux) {
            Copy-Item (Join-Path $ZLIBInstallDir 'lib' 'libz.so')      -Destination (Join-Path $NetCDFBinDir 'libz.so')
        }
        else {
            Copy-Item (Join-Path $ZLIBInstallDir 'lib' 'libz.dylib')      -Destination (Join-Path $NetCDFBinDir 'libz.dylib')
        }
    }
    elseif ($IsWindows) {
        Copy-Item (Join-Path $ZLIBInstallDir 'lib' 'zlibstatic.lib') -Destination (Join-Path $NetCDFLibDir 'zlib.lib')
        Copy-Item (Join-Path $HDF5InstallDir 'lib' 'libhdf5.lib')    -Destination (Join-Path $NetCDFLibDir 'hdf5-static.lib')
        Copy-Item (Join-Path $HDF5InstallDir 'lib' 'libhdf5_hl.lib') -Destination (Join-Path $NetCDFLibDir 'hdf5_hl-static.lib')
    }

    # Back-up  $NetCDFInstallDir/lib/cmake/netCDF/netCDFTargets.cmake before modifying it.
    # For reference, the back-up will not be deleted.
    $NeCDFCMakeTargets = (Join-Path $NetCDFLibDir 'cmake' 'netCDF' 'netCDFTargets.cmake')
    Copy-Item $NeCDFCMakeTargets -Destination (Join-Path $NetCDFLibDir 'cmake' 'netCDF' 'netCDFTargets.cmake.original')

    # Find the line in set_target_properties where INTERFACE_LINK_LIBRARIES is set
    $Content = (Get-Content -Path $NeCDFCMakeTargets)
    $Line = $Content | Select-String "INTERFACE_LINK_LIBRARIES" | Select-Object -ExpandProperty Line
    if ($null -eq $Line ) {
        Write-Error ('Cannot modify { 0 }. Please contact the developers.' -f $NeCDFCMakeTargets)
        Invoke-Terminate
        Exit
    }

    # Replace the line above to list the public interface libararies in ${_IMPORT_PREFIX}/lib (libs copied and renamed above)
    if ($IsLinux) {
        $NewLine = '  INTERFACE_LINK_LIBRARIES "dl;${_IMPORT_PREFIX}/lib/libhdf5_hl.a;${_IMPORT_PREFIX}/lib/libhdf5.a;${_IMPORT_PREFIX}/lib/libz.a;${_IMPORT_PREFIX}/bin/libz.so"'
    }
    elseif ($IsMacOS) {
        $NewLine = $Line
        $NewLine.Replace($HDF5InstallDir, ${_IMPORT_PREFIX})
        $NewLine.Replace((Join-Path $ZLIBInstallDir "lib"), (Join-Path ${_IMPORT_PREFIX} "bin"))
    }
    elseif ($IsWindows) {
        $NewLine = '  INTERFACE_LINK_LIBRARIES "${_IMPORT_PREFIX}/lib/hdf5_hl-static.lib;${_IMPORT_PREFIX}/lib/hdf5-static.lib;${_IMPORT_PREFIX}/lib/zlib.lib"'
    }
    $Content.Replace($Line, $NewLine) | Set-Content $NeCDFCMakeTargets
}

Invoke-Post-Build-Steps

Write-Host ( `
        "`nCongratulations! You did it!`n" `
        + "`nYou must now either:`n" `
        + (" - add {0} to your path before building MeshKernel, or`n" -f $NetCDFInstallDir) `
        + (" - configure the MeshKernel build with -DCMAKE_PREFIX_PATH={0}." -f $NetCDFInstallDir) `
)

Invoke-Terminate

if ($Clean) {
    Remove-Item -Recurse -Force $WorkDir
}
