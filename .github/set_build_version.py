#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path

VERSION_LINE_RE = re.compile(r'^(version\s*=\s*")([^"]+)("\s*)$', re.MULTILINE)


def stamp_version(version: str, build_number: str) -> str:
    public = version.split("+", maxsplit=1)[0]
    return f"{public}+build.{build_number}"


def update_pyproject(pyproject_path: Path, build_number: str) -> str:
    content = pyproject_path.read_text(encoding="utf-8")
    match = VERSION_LINE_RE.search(content)
    if match is None:
        raise RuntimeError("Could not find [project].version in pyproject.toml")

    original_version = match.group(2)
    new_version = stamp_version(original_version, build_number)
    updated = VERSION_LINE_RE.sub(
        lambda m: f"{m.group(1)}{new_version}{m.group(3)}",
        content,
        count=1,
    )
    pyproject_path.write_text(updated, encoding="utf-8")
    return new_version


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stamp pyproject version with a CI build number"
    )
    parser.add_argument(
        "--pyproject",
        type=Path,
        required=True,
        help="Path to pyproject.toml",
    )
    parser.add_argument(
        "--build-number",
        type=str,
        required=True,
        help="Build number for the +build.N local version suffix",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    new_version = update_pyproject(args.pyproject, args.build_number)
    print(f"Stamped build version: {new_version}")


if __name__ == "__main__":
    main()
