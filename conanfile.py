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
