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


def build_index(wheels_dir: Path, output_dir: Path) -> None:
    package_dir = output_dir / "simple" / "prexsyn-engine"
    package_dir.mkdir(parents=True, exist_ok=True)

    wheel_files = sorted(wheels_dir.glob("*.whl"))
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
    write_page(package_dir / "index.html", links or "    <!-- no wheels found -->")


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
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_index(args.wheels, args.output)


if __name__ == "__main__":
    main()
