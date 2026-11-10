import os

from conan import ConanFile
from conan.tools.env import Environment
from conan.tools.cmake import CMakeDeps, CMakeToolchain


class MeshKernelConan(ConanFile):
    name = "meshkernel"
    version = "1.0"

    settings = "os", "arch", "compiler", "build_type"

    requires = "boost/1.78.0"

    default_options = {
        "boost/*:header_only": True,
    }

    def generate(self):
        # Generate CMake dependency configuration.
        deps = CMakeDeps(self)
        deps.generate()

        # Generate CMake toolchain.
        tc = CMakeToolchain(self)
        tc.generate()

        # ------------------------------------------
        # Linux: GCC
        # ------------------------------------------
        if self.settings.os == "Linux":
            gcc_path = "/opt/apps/gcc/12.2.0_gcc12.2.0/bin"

            env = Environment()
            env.define("CC", os.path.join(gcc_path, "gcc"))
            env.define("CXX", os.path.join(gcc_path, "g++"))

            env.vars(self).save_script("gcc12")

        # ------------------------------------------
        # Windows: Microsoft Visual C++
        # ------------------------------------------
        elif self.settings.os == "Windows":
            # MSVC is initialized by the Visual Studio
            # CMake generator.
            #
            # Generator and compiler version are defined
            # in the Conan profile.
            pass

        # ------------------------------------------
        # macOS: Apple Clang
        # ------------------------------------------
        elif self.settings.os == "Macos":
            env = Environment()
            env.define("CC", "/usr/bin/clang")
            env.define("CXX", "/usr/bin/clang++")

            env.vars(self).save_script("appleclang")

