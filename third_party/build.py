#!/usr/bin/env python3
"""Build the static third-party libraries used by prexsyn-engine.

This replaces the platform-specific shell bootstrap during CMake builds.  It is
deliberately dependency-free so it can run in isolated wheel build environments.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys


BOOST_VERSION = "boost-1.86.0"
RDKIT_VERSION = "Release_2025_09_5"
BOOST_LIBRARIES = (
    "serialization",
    "iostreams",
    "program_options",
    "json",
    "stacktrace",
    "format",
    "graph",
    "rational",
    "flyweight",
    "math",
    "property_tree",
    "crc",
    "multiprecision",
    "assign",
    "multi_array",
)


def run(command: list[str], *, cwd: Path) -> None:
    print("+", subprocess.list2cmdline(command), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def clone(url: str, branch: str, destination: Path) -> None:
    if destination.is_dir():
        return
    run(
        ["git", "clone", url, "--branch", branch, "--depth", "1", str(destination)],
        cwd=destination.parent,
    )


def build_boost(source_root: Path, prefix: Path) -> None:
    marker = prefix / f".prexsyn-{BOOST_VERSION}"
    if marker.exists():
        print(f"Boost {BOOST_VERSION} is already installed in {prefix}")
        return

    source = source_root / "boost"
    clone("https://github.com/boostorg/boost.git", BOOST_VERSION, source)
    run(
        ["git", "submodule", "update", "--depth", "1", "--init", "tools/boostdep"],
        cwd=source,
    )
    for library in BOOST_LIBRARIES:
        run(
            [
                "git",
                "submodule",
                "update",
                "--depth",
                "1",
                "--init",
                f"libs/{library}",
            ],
            cwd=source,
        )
        run(
            [
                sys.executable,
                "tools/boostdep/depinst/depinst.py",
                "-X",
                "test",
                "-g",
                "--depth 1",
                library,
            ],
            cwd=source,
        )

    prefix_arg = f"--prefix={prefix}"
    if os.name == "nt":
        # Boost 1.86 detects Visual Studio 18 as the unsupported "vcunk"
        # toolset. CI uses Visual Studio 2022, so select its toolset explicitly.
        run(["cmd.exe", "/d", "/c", "bootstrap.bat", "vc143"], cwd=source)
        b2 = source / "b2.exe"
    else:
        bootstrap = ["./bootstrap.sh", prefix_arg]
        toolset = os.environ.get("BOOST_TOOLSET")
        if toolset:
            bootstrap.append(f"--with-toolset={toolset}")
        run(bootstrap, cwd=source)
        b2 = source / "b2"

    b2_command = [
        str(b2),
        "variant=release",
        "link=static",
        "threading=multi",
        prefix_arg,
    ]
    if os.name == "nt":
        b2_command.append("runtime-link=shared")
    else:
        b2_command.extend(("cxxflags=-fPIC", "cflags=-fPIC"))
    b2_command.append("install")
    run(b2_command, cwd=source)
    marker.touch()


def build_rdkit(
    source_root: Path,
    prefix: Path,
    generator: str,
    platform: str | None,
    toolset: str | None,
    additional_prefixes: list[str],
    static_boost: bool,
) -> None:
    marker = prefix / f".prexsyn-rdkit-{RDKIT_VERSION}"
    if marker.exists():
        print(f"RDKit {RDKIT_VERSION} is already installed in {prefix}")
        return

    source = source_root / "rdkit"
    build = source / "build-prexsyn"
    clone("https://github.com/rdkit/rdkit.git", RDKIT_VERSION, source)

    prefix_path = ";".join((str(prefix), *additional_prefixes))
    configure = [
        "cmake",
        "-S",
        str(source),
        "-B",
        str(build),
        "-G",
        generator,
        "-DCMAKE_BUILD_TYPE=Release",
        f"-DCMAKE_PREFIX_PATH={prefix_path}",
        f"-DCMAKE_INSTALL_PREFIX={prefix}",
        "-DRDK_INSTALL_INTREE=OFF",
        "-DRDK_INSTALL_STATIC_LIBS=ON",
        "-DRDK_BUILD_CPP_TESTS=OFF",
        "-DRDK_BUILD_PYTHON_WRAPPERS=OFF",
        "-DRDK_BUILD_CHEMDRAW_SUPPORT=OFF",
        "-DRDK_BUILD_CAIRO_SUPPORT=OFF",
        "-DRDK_BUILD_FREETYPE_SUPPORT=OFF",
        "-DRDK_INSTALL_COMIC_FONTS=OFF",
        "-DRDK_BUILD_DESCRIPTORS3D=OFF",
    ]
    if static_boost:
        configure.append("-DBoost_USE_STATIC_LIBS=ON")
    if platform:
        configure.extend(("-A", platform))
    if toolset:
        configure.extend(("-T", toolset))
    run(configure, cwd=source_root)
    run(
        [
            "cmake",
            "--build",
            str(build),
            "--config",
            "Release",
            "--parallel",
            str(os.cpu_count() or 2),
        ],
        cwd=source_root,
    )
    run(
        ["cmake", "--install", str(build), "--config", "Release"],
        cwd=source_root,
    )
    marker.touch()


def clean(prefix: Path) -> None:
    """Remove only artifacts created by this script beneath *prefix*."""
    for directory in ("bin", "include", "lib", "share"):
        path = prefix / directory
        if path.is_dir():
            print(f"Removing {path}")
            shutil.rmtree(path)
    for marker in prefix.glob(".prexsyn-*"):
        if marker.is_file():
            print(f"Removing {marker}")
            marker.unlink()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", type=Path)
    parser.add_argument("--generator")
    parser.add_argument("--platform")
    parser.add_argument("--toolset")
    parser.add_argument("--cmake-prefix", action="append", default=[])
    parser.add_argument(
        "--rdkit-static-boost",
        action="store_true",
        help="Link RDKit against static Boost libraries",
    )
    parser.add_argument("--skip-boost", action="store_true")
    parser.add_argument("--skip-rdkit", action="store_true")
    parser.add_argument("--clean", action="store_true")
    args = parser.parse_args()
    if not args.clean and not args.generator:
        parser.error("--generator is required unless --clean is used")
    return args


def main() -> None:
    args = parse_args()
    prefix = args.prefix.resolve()
    prefix.mkdir(parents=True, exist_ok=True)
    if args.clean:
        clean(prefix)
        return

    source_root = prefix / "src"
    source_root.mkdir(exist_ok=True)

    if not args.skip_boost:
        build_boost(source_root, prefix)
    if not args.skip_rdkit:
        build_rdkit(
            source_root,
            prefix,
            args.generator,
            args.platform,
            args.toolset,
            args.cmake_prefix,
            args.rdkit_static_boost,
        )


if __name__ == "__main__":
    main()
