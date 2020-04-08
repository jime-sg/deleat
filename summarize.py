#!/usr/bin/env python3
"""
# TODO
@author: Jimena Solana
"""

from argparse import ArgumentParser

import circplot


if __name__ == "__main__":
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="summarize",
        description=""  # FIXME
    )

    circplot.plot(
        genbank="/home/jimena/Bartonella/NC_005955_wregions.gb",
        out_file="/home/jimena/Dropbox/TFM/paso_circos/test4",
        out_fmt="png"
    )
