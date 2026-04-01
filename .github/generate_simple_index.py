#!/usr/bin/env python3

from __future__ import annotations

import argparse
from html import escape
from pathlib import Path


def write_page(path: Path, body: str) -> None:
    path.write_text(
        f"<!doctype html>\n<html>\n  <body>\n{body}\n  </body>\n</html>\n",
        encoding="utf-8",
    )


def extract_versions(wheel_files: list[Path]) -> list[str]:
    versions: list[str] = []
    for wheel in wheel_files:
        parts = wheel.name.split("-")
        if len(parts) < 2:
            continue
        version = parts[1]
        if version not in versions:
            versions.append(version)
    return versions


def build_index(wheels_dir: Path, output_dir: Path, build_number: str) -> None:
    package_dir = output_dir / "simple" / "prexsyn-engine"
    package_dir.mkdir(parents=True, exist_ok=True)

    wheel_files = sorted(wheels_dir.glob("*.whl"))
    versions = extract_versions(wheel_files)
    if not versions:
        version_label = "unknown"
    elif len(versions) == 1:
        version_label = versions[0]
    else:
        version_label = ", ".join(versions)
    version_tag = escape(f"{version_label}+build.{build_number}")

    links = "\n".join(
        f'    <a href="{escape(wheel.name)}">{escape(wheel.name)}</a><br/>'
        for wheel in wheel_files
    )

    write_page(
        output_dir / "index.html",
        '    <a href="simple/">simple/</a>',
    )
    write_page(
        output_dir / "simple" / "index.html",
        '    <a href="prexsyn-engine/">prexsyn-engine/</a>',
    )
    package_body = f"    <p>version tag: {version_tag}</p>\n" + (
        links or "    <!-- no wheels found -->"
    )
    write_page(package_dir / "index.html", package_body)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate a simple PyPI index")
    parser.add_argument(
        "--wheels",
        type=Path,
        required=True,
        help="Directory containing wheel files",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Directory where the index should be written",
    )
    parser.add_argument(
        "--build-number",
        type=str,
        required=True,
        help="Build number used in the rendered version tag",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_index(args.wheels, args.output, args.build_number)


if __name__ == "__main__":
    main()
