import os
import pathlib
import sys

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name: str, sourcedir: str = "") -> None:
        super().__init__(name, sources=[])
        self.sourcedir = os.fspath(pathlib.Path(sourcedir).resolve())


class CMakeBuild(build_ext):
    def run(self) -> None:
        for ext in self.extensions:
            self.build_cmake(ext)

    def build_cmake(self, ext: CMakeExtension) -> None:
        cwd = pathlib.Path().absolute()

        build_temp = pathlib.Path(self.build_temp)
        print(f"build_temp: {build_temp}")
        build_temp.mkdir(parents=True, exist_ok=True)
        extdir = pathlib.Path(self.get_ext_fullpath(ext.name)).parent.resolve()

        cmake_args = [
            "-DCMAKE_PREFIX_PATH=" + sys.prefix,
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + str(extdir),
            "-DCMAKE_BUILD_TYPE=Release",
            "-G",
            "Ninja",
        ]

        build_args = ["--config", "Release", "--", "-j4"]

        os.chdir(str(build_temp))
        self.spawn(["cmake", ext.sourcedir] + cmake_args)
        if not self.dry_run:
            self.spawn(["cmake", "--build", "."] + build_args)
        os.chdir(str(cwd))

        # Find all submodules. Move submodules into subpackage directory
        # For example: a__b.cpython-313-x86_64-linux-gnu.so -> a/b.cpython-313-x86_64-linux-gnu.so
        for file in extdir.glob("**/*.so"):
            parts = file.name.split("__")
            if len(parts) > 1:
                subpackage_dir = extdir.joinpath(*parts[:-1])
                subpackage_dir.mkdir(parents=True, exist_ok=True)
                file.rename(subpackage_dir.joinpath(parts[-1]))


setup(
    packages=find_packages(include=["prexsyn_engine", "prexsyn_engine.*"]),
    ext_modules=[CMakeExtension("prexsyn_engine.*")],
    cmdclass={"build_ext": CMakeBuild},
)
