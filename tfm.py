#!/usr/bin/env python3
"""tfm.py
@author: Jimena Solana
"""

from sys import argv
import runpy


STEPS = [
    "nonessential-genes",
    "define-deletions",
    "revise-deletions",
    "design-all-primers",
    "summarize"
]


def print_usage():
    print("usage: tfm <step name> <step arguments>")
    print("valid steps are:")
    for i, s in enumerate(STEPS):
        print("\t%d. %s" % (i + 1, s))
    raise SystemExit


if __name__ == "__main__":
    if len(argv) < 2:
        print_usage()
    step = argv[1]
    if step not in STEPS:
        if step not in ("-h", "--help"):
            print("error: unrecognized option '%s'" % step)
        print_usage()

    else:
        runpy.run_module(step.replace("-", "_"), run_name="__main__")
