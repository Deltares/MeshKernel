# MeshKernel
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=Deltares_Grid_Editor_back-end&metric=alert_status)](https://sonarcloud.io/dashboard?id=Deltares_Grid_Editor_back-end)

MeshKernel is a library for creating and editing meshes.
It supports 1D & 2D unstructured meshes as well as curvilinear meshes.

Wrappers for the following languages are available:
- [Python](https://github.com/Deltares/MeshKernelPy)
- [.NET](https://github.com/Deltares/MeshKernelNET)

The library is separated in an API namespace (MeshKernelApi) used for communication with the client and a backend namespace (MeshKernel) where the algorithms are implemented. 
The API namespace contains several structures used as parameters for the API methods (see API usage section). 
These structures must be mirrored in the client application and filled with appropriate values.

## Examples

1. Creating a triangular mesh inside a polygon

In this example a mesh is created by discretizing the polygon perimeter with the desired edge length.

![alt tag](docs/images/TriangularMeshInPolygon.jpg)

2. Mesh orthogonalization

Finite volume staggered flow solvers require the mesh to be as much orthogonal as possible. 
MeshKernel provides an algorithm to adapt the mesh and achieve a good balance between mesh orthogonality and smoothness.

![alt tag](docs/images/MeshOrthogonalization.jpg)

3. Curvilinear mesh generation

Curvilinear meshes for rivers can be generated using splines.

![alt tag](docs/images/OrthogonalCurvilinearGrid.jpg)

4. Mesh refinement

A mesh can be refined in areas based on samples or polygon selections. 

![alt tag](docs/images/GridRefinement.jpg)


## Shared library dependencies (Linux)
- libgomp

## Build Requirements

The requirements are:
- Git
- CMake 3.23 or higher
- A C++20 compatible compiler
- The Boost libraries
- NetCDF
- Doxygen (optional)

## Dependencies

### Boost
- Under Windows, precompiled boost binaries (with MSVC compiler) can be downloaded [here](https://sourceforge.net/projects/boost/files/boost-binaries/). Alternatively, the source code is available [here](https://sourceforge.net/projects/boost/files/boost/) and the installation instructions can be found [here](https://www.boost.org/doc/libs/1_74_0/more/getting_started/windows.html).
- Under Linux, Boost can be either obtained from the package repository of the used Linux distribution or built from [source](https://sourceforge.net/projects/boost/files/boost/) following these [instructions](https://www.boost.org/doc/libs/1_74_0/more/getting_started/unix-variants.html).

 Once installed, add the boost environmental variables accordingly. For example, if version 1.74.0 in installed in `C:\Apps`, set the environment variable:
  ```
  BOOST_INCLUDEDIR=C:\Apps\boost_1_74_0
  ```

### NetCDF
The NetCDF static library is required for building MeshKernel. Optionally, a PowerShell [script](scripts/install_netcdf_static.ps1) is made available for cloning and building [NetCDF](https://github.com/Unidata/netcdf-c) and its dependencies ([HDF5](https://github.com/HDFGroup/hdf5), [ZLIB](https://github.com/madler/zlib), [Curl](https://github.com/curl/curl), and [m4](https://sourceforge.net/projects/gnuwin32/files/m4/)). It can be used under both Windows and Linux.

To run the script in a PowerShell session, use
```powershell
.\install_netcdf_static.ps1 -WorkDir '/path/to/work/directory' -InstallDir '/path/to/install/directory' -BuildType 'Release' -ParallelJobs 10 -GitTags @{zlib = 'v1.2.13'; curl = 'curl-7_88_1';  hdf5 = 'hdf5-1_14_0';  netcdf_c = 'v4.9.1'}
```

with `/path/to/work/directory` and `/path/to/install/directory` replaced with valid paths.

For more information regarding the script's options above, use
```powershell
Get-Help .\install_netcdf_static.ps1 -Detailed
```

Upon successful installation, to build MeshKernel successfully, it is important to either
- add the path to the install directory to the system path, or
- configure the MeshKernel build with `-DCMAKE_PREFIX_PATH=/path/to/install/directory`.

**Note:** Additional dependencies may be required depending on the system configuration:
- Windows: [Perl](https://strawberryperl.com/)
- Linux: m4, OpenSSL, Curl, and [PowerShell](https://learn.microsoft.com/en-us/powershell/scripting/install/installing-powershell-on-linux?view=powershell-7.3). All dependencies can be installed from the repository of the used Linux distribution.
  
## Configuring and Building MeshKernel
Follow the steps below to configure, build and install MeshKernel.
### Steps
1.  To configure under Windows with Visual Studio, a solution is generated using

    ```powershell
    cmake -S <path-to-source-dir> -B <path-to-build-dir> -G "Visual Studio 16 2019" --install-prefix <path-to-install-dir>
    ```
    This example uses Visual Studio 16 2019. A different version can be specified. In the above
    - `<path-to-source-dir>` is the path to the MeshKernel source tree (the top-level directory containing source files provided by the project).
    - `<path-to-build-dir>` is the path to the directory where MeshKernel is to be built.
    - `<path-to-install-dir>` is the path top the directory where MeshKernel is to be installed.

    Under Linux, the generator option is omitted:
     ```powershell
    cmake -S <path-to-source-dir> -B <path-to-build-dir> --install-prefix <path-to-install-dir>
    ```

2.  To build the project's targets, use:
    ```powershell
    cmake --build <path-to-build-dir> --config <cfg> --parallel <jobs>
    ```
    where
    - `<cfg>` is the build type (`Debug`, `Release`, `RelWithDebInfo` and `MinSizeRel`), see `CMAKE_BUILD_TYPE`.
    - `<jobs>` is the maximum number of concurrent processes to use when building.

    Note: this builds the documentation by default in `<path-to-build-dir>/docs/html`.

3.  To install, use:
    ```powershell
    cmake --install <path-to-build-dir>
    ```
    or 
    ```powershell
    cmake --install <path-to-build-dir> --prefix [<path-to-install-dir>]
    ```
    to override the installation path specified during configuration (step 1).


### Additional configuration options
MeshKernel can be configured with a set of options, which are summarized in the table below.

| Option | Description | Type | Default value | Notes |
|------- |-------------|------|---------------|-------|
| ENABLE_UNIT_TESTING | Enables building the unit test executables | Bool | ON  | |
| ENABLE_BENCHMARKING | Enables building the benchmark executable | Bool | OFF | |
| ENABLE_BENCHMARKING_MEM_REPORT | Enables reporting memory usage statistics | Bool | OFF | Applicable only when ENABLE_BENCHMARKING is ON, ignored otherwise |
| ENABLE_CODE_COVERAGE | Generates code coverage report | Bool | OFF | Valid only under Linux when a GNU compiler is used (requires gcov), ignored otherwise |








