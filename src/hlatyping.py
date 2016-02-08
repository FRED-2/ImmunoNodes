#!/usr/bin/env python
"""
Commandline tool for hla typing prediction

usage: hlatyping.py [-h] [-m {optitype,athlates,seq2hla,polysolver}]
                    [-v VERSION] -i INPUT [-p PAIRED] [-r {rna,dna}] -o OUTPUT

Commandline tool for HLA typing

optional arguments:
  -h, --help            show this help message and exit
  -m, --method {optitype,athlates,seq2hla,polysolver}
                        The name of the prediction method
  -v VERSION, --version VERSION
                        The version of the prediction method
  -i INPUT, --input INPUT
                        Path to the input file
  -p PAIRED, --paired PAIRED
                        Additional input for paired-end typing
  -r {rna,dna}, --reference {rna,dna}
                        The reference type to use
  -o OUTPUT, --output OUTPUT
                        Path to the output file

"""
import sys
import tempfile
import argparse

from Fred2.HLAtyping import HLATypingFactory


def main():
#Specify CTD interface
    # Every CTD Model has to have at least a name and a version, plus any of the optional attributes below them.
    model = argparse.ArgumentParser(
        description='Commandline tool for HLA typing',
        )

    model.add_argument('-m',
        '--method',
        type=str,
        choices=HLATypingFactory.available_methods().keys(),
        default="optitype",
        help='The name of the prediction method'
        )

    model.add_argument('-v',
        '--version',
        type=str,
        default="",
        help='The version of the prediction method'
        )

    model.add_argument('-i',
        '--input',
        type=str,
        required=True,
        help='Path to the input file'
        )

    model.add_argument('-p',
        '--paired',
        type=str,
        default="",
        help="Additional input for paired-end typing"
    )

    model.add_argument('-r',
        '--reference',
        type=str,
        choices=["rna", "dna"],
        default="dna",
        help='The reference type to use'
        )

    model.add_argument('-o',
        '--output',
        type=str,
        required=True,
        help='Path to the output file'
        )

    args = model.parse_args()

    version = "" if args.version == "" else args.version
    options = "--"+args.reference if args.paired == "" else args.paired+" "+"--"+args.reference

    genotype = HLATypingFactory(args.method).predict(args.input, "/tmp/", options=options)
    with open(args.output, "w") as f:
        f.write("\n".join("HLA-"+a.name for a in genotype))

    return 0


if __name__ == "__main__":
    sys.exit(main())