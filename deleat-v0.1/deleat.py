#!/usr/bin/env python3
"""deleat.py

    DELEAT: deletion design by essentiality analysis tool
    Version 0.1

@author: Jimena Solana
"""

from sys import argv
import runpy


STEPS = [
    "predict-essentiality",
    "define-deletions",
    "revise-deletions",
    "summarise",
    "design-all-primers",
    "design-primers"
]


def print_usage():
    print("usage: deleat <step name> <step arguments>")
    print("valid steps are:")
    print(
        "  1. predict-essentiality"
        "\n  2. define-deletions"
        "\n  3. revise-deletions"
        "\n  4. summarise"
        "\n  5. design-all-primers / design-primers"
    )
    raise SystemExit


if __name__ == "__main__":
    if len(argv) < 2:
        print_usage()
    step = argv[1]
    if step not in STEPS:
        if step not in ("-h", "--help"):
            print("error: unrecognised option '%s'" % step)
        print_usage()

    else:
        runpy.run_module(step.replace("-", "_"), run_name="__main__")
