# MeshKernel: Conan integration for Boost

This document describes the Conan 2 integration used by MeshKernel to
obtain **Boost 1.78.0** from Conan Center. Other dependencies, notably
the compiler, CMake and netCDF, remain provided by the existing platform
toolchains and environment modules.

The examples use a Release build. Adjust profile names and toolchain
paths to match the machine being used.

## 1. Repository layout

Keep the Conan recipe in the repository root. The provider and platform
profiles live under `conan/`.

``` text
MeshKernel/
├── CMakeLists.txt
├── conanfile.py
└── conan/
    ├── cmake/
    │   └── conan_provider.cmake
    └── profiles/
        ├── gcc12-release
        ├── macos-clang-release
        └── windows-msvc-2022
```

## 2. Enable the Conan dependency provider

The Conan dependency provider must be loaded before the first
`project()` call. Add the following near the beginning of the root
`CMakeLists.txt`:

``` cmake
cmake_minimum_required(VERSION 3.24)

set(CMAKE_PROJECT_TOP_LEVEL_INCLUDES
    "${CMAKE_CURRENT_SOURCE_DIR}/conan/cmake/conan_provider.cmake")

set(MESHKERNEL_VERSION 8.3.0)

project(
    MeshKernel
    VERSION ${MESHKERNEL_VERSION}
    # Keep the remaining project arguments here.
)
```

CMake 3.24 or newer is required for this dependency-provider
integration.

### Obtain the provider

The repository includes a pinned copy of `conan_provider.cmake` from the
`cmake-conan` 0.19.0 release.

**Linux**

``` bash
mkdir -p conan/cmake
curl -L \
    https://raw.githubusercontent.com/conan-io/cmake-conan/0.19.0/conan_provider.cmake \
    -o conan/cmake/conan_provider.cmake
```

**Windows PowerShell**

``` powershell
New-Item -ItemType Directory -Force .\conan\cmake | Out-Null

Invoke-WebRequest -UseBasicParsing `
    -Uri "https://raw.githubusercontent.com/conan-io/cmake-conan/0.19.0/conan_provider.cmake" `
    -OutFile .\conan\cmake\conan_provider.cmake
```

Commit the provider file to the repository so that all developers and
build agents use the same pinned version.

## 3. Conan recipe

The root-level `conanfile.py` declares Boost as a header-only
dependency:

``` python
from conan import ConanFile


class MeshKernelConan(ConanFile):
    name = "meshkernel"
    version = "1.0"

    settings = "os", "arch", "compiler", "build_type"

    requires = "boost/1.78.0"

    default_options = {
        "boost/*:header_only": True,
    }

    generators = "CMakeConfigDeps"
```

`CMakeConfigDeps` generates the CMake package configuration files
consumed by `find_package()`. When using the CMake dependency provider,
use this generator rather than the older explicit `CMakeDeps` plus
`CMakeToolchain` workflow.

Do not add a `generate()` method that calls `CMakeDeps` or
`CMakeToolchain` when using the recipe above. The declared generator is
run automatically by Conan.

## 4. Conan profiles

Profiles describe the target platform and compiler. Use the profile
appropriate to the machine being built on.

### Linux: GCC 12

`conan/profiles/gcc12-release`:

``` ini
[settings]
arch=x86_64
build_type=Release
compiler=gcc
compiler.cppstd=gnu17
compiler.libcxx=libstdc++11
compiler.version=12
os=Linux
```

### macOS: Apple Clang

`conan/profiles/macos-clang-release`:

``` ini
[settings]
arch=armv8
build_type=Release
compiler=apple-clang
compiler.cppstd=gnu17
compiler.libcxx=libc++
compiler.version=<detected-version>
os=Macos
```

The compiler version is machine-dependent and must be replaced with the
version detected on the target Mac. This profile is a starting point and
has not been validated in the source notes.

### Windows: Visual Studio 2022 or MSBuild Tools 2022

`conan/profiles/windows-msvc-2022`:

``` ini
[settings]
arch=x86_64
build_type=Release
compiler=msvc
compiler.cppstd=20
compiler.runtime=dynamic
compiler.runtime_type=Release
compiler.version=194
os=Windows

[conf]
tools.cmake.cmaketoolchain:generator=Visual Studio 17 2022
```

The profile assumes MSVC 19.4x / Visual Studio 2022. Confirm that the
selected compiler version and runtime match the active Visual Studio
toolchain.

## 5. General Conan commands

The commands below are shown for Linux. The Conan commands are the same
on Windows, apart from shell syntax and installation paths.

### Install and check Conan

``` bash
python -m pip install --user conan
conan --version
```

Ensure that the directory containing the Conan executable is on `PATH`.
On Windows, update `$env:Path` as appropriate.

The examples assume Conan 2 is installed and available as `conan`.

### Configure Conan on a new machine

On a machine without the Deltares profiles, detect a default profile:

``` bash
conan profile detect --force
```

This is normally unnecessary on Deltares build machines where the
repository profiles are used.

### Inspect remotes and available Boost versions

``` bash
conan remote list
conan list "boost/*" -r=conancenter
```

### Clear the Conan cache

To remove all cached recipes and packages:

``` bash
conan remove "*" --confirm
```

This is destructive to the local Conan cache, although it does not
remove the profiles and remote configuration under the Conan home
directory. Use it only when a complete cache reset is intended. It is
not required for normal builds.

To remove only Boost:

``` bash
conan remove "boost/*" --confirm
```

## 6. Linux build

The Linux build uses Environment Modules for GCC, CMake and netCDF.
Conan supplies Boost.

### Load the toolchain

``` bash
module purge
module load gcc/12.2.0_gcc12.2.0
module load cmake/3.26.4_gcc12.2.0
module load netcdf/4.9.2_gcc12.2.0

export CC=/opt/apps/gcc/12.2.0_gcc12.2.0/bin/gcc
export CXX=/opt/apps/gcc/12.2.0_gcc12.2.0/bin/g++
```

Setting `CC` and `CXX` explicitly ensures CMake uses the intended GCC
installation rather than a system compiler.

### Configure Release

``` bash
cmake -S . -B build-conan-release \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PROJECT_TOP_LEVEL_INCLUDES="$PWD/conan/cmake/conan_provider.cmake" \
    -DCMAKE_FIND_PACKAGE_PREFER_CONFIG=ON \
    -DCMAKE_PREFIX_PATH="/opt/apps/netcdf/4.9.2_gcc12.2.0" \
    -DCMAKE_INSTALL_PREFIX="$PWD/install-conan-release" \
    -DCONAN_HOST_PROFILE="$PWD/conan/profiles/gcc12-release" \
    -DCONAN_BUILD_PROFILE="$PWD/conan/profiles/gcc12-release" \
    -DENABLE_UNIT_TESTING=OFF \
    -DENABLE_BENCHMARKING=OFF
```

The provider invokes Conan during CMake configuration when the project
first requests a package with `find_package()`. Conan resolves Boost
using the selected profile and writes generated CMake package files
under the build directory.

### Build

Load the environment script generated by the Conan recipe only if the
recipe defines one. For the minimal `CMakeConfigDeps` recipe above, no
`gcc12.sh` environment script is generated; the compiler is selected
before configuration through `CC`, `CXX` and the loaded module
environment.

``` bash
cmake --build build-conan-release --config Release --parallel 8
```

### Install

``` bash
cmake --install build-conan-release
```

The install prefix is `install-conan-release` as specified during
configuration.

## 7. Windows build

The Windows build uses Visual Studio 2022 and the installed netCDF
distribution. Conan supplies Boost.

### Select the profiles

In PowerShell:

``` powershell
$env:CONAN_HOST_PROFILE = "$PWD/conan/profiles/windows-msvc-2022"
$env:CONAN_BUILD_PROFILE = "$PWD/conan/profiles/windows-msvc-2022"
```

### Configure Release

``` powershell
cmake -S . -B build-conan-release `
    -G "Visual Studio 17 2022" `
    -A x64 `
    --install-prefix "D:/development/MeshKernel/install-conan-release" `
    -DCMAKE_PROJECT_TOP_LEVEL_INCLUDES="$PWD/conan/cmake/conan_provider.cmake" `
    -DCMAKE_FIND_PACKAGE_PREFER_CONFIG=ON `
    "-DCMAKE_PREFIX_PATH=C:/Program Files (x86)/netCDF" `
    -DENABLE_UNIT_TESTING=OFF `
    -DENABLE_BENCHMARKING=OFF
```

Visual Studio is a multi-configuration generator, so the configuration
is selected when building rather than with `CMAKE_BUILD_TYPE`.

### Build

``` powershell
cmake --build build-conan-release --config Release --parallel 8
```

### Install

``` powershell
cmake --install build-conan-release --config Release
```

## 8. Notes and troubleshooting

-   **Keep the Conan recipe in the repository root.** The unmodified
    `conan_provider.cmake` invokes Conan with the CMake source directory
    as the recipe path. It does not automatically search a `conan/`
    subdirectory for `conanfile.py`.
-   **Use `CMakeConfigDeps` with the dependency provider.** Do not mix
    this with the older explicit `CMakeDeps`/`CMakeToolchain` generation
    workflow.
-   **Use a clean build directory when changing provider or toolchain
    settings.** CMake caches toolchain and configuration values. For
    example, remove `build-conan-release` and configure it again.
-   **Keep netCDF outside Conan for this integration.** The current
    scope is Boost only; netCDF and the compiler toolchain continue to
    come from the machine environment.
-   **The macOS profile is provisional.** Its compiler version and
    settings need validation on the intended Mac environment.
-   **The provider version is pinned.** If it is upgraded, review the
    change and test the integration on each supported platform.
