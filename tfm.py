#!/usr/bin/env python3
"""
@author: Jimena Solana
"""

from argparse import ArgumentParser
import runpy


if __name__ == "__main__":
    parser = ArgumentParser(
        description=""  # FIXME
    )
    parser.add_argument(
        dest="step",
        help="",  # FIXME
        choices=(
            "nonessential-genes",
            "nonessential-regions",
            "deletion-coords",
            "design-primers",
            "summarize"
        )
    )
    args = parser.parse_args()
    step = args.step

    if step == "nonessential-genes":
        runpy.run_module("nonessential_genes")  # FIXME
    elif step == "nonessential-regions":
        runpy.run_module("nonessential_regions")  # FIXME
    elif step == "deletion-coords":
        runpy.run_module("integrate_scores")  # FIXME
    elif step == "design-primers":
        runpy.run_module("primer_design")  # FIXME
    elif step == "summarize":
        runpy.run_module("summarize")  # FIXME


"""
a.py:
#!/usr/bin/env python3
import runpy
runpy.run_module("b", run_name="__main__")

b.py:
#!/usr/bin/env python3
from sys import argv
if __name__ == "__main__":
    print("lo que me ha llegado es:", argv[1])
"""
