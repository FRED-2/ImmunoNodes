#!/usr/bin/env python
"""
Command line tool for generation of trie data structure for a set of proteins

usage: distance2selfgeneration.py [-h] -i INPUT [-l LENGTH] -o OUTPUT
                                  [-b BLOSUM]

Generation of trie data structure for a set of proteins for distance to self
calculation

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Protein sequences in fasta file format. All peptides
                        of the specified length will be generated for the trie
                        generation.
  -l LENGTH, --length LENGTH
                        Specifies the length of peptides which will be
                        generated.
  -o OUTPUT, --output OUTPUT
                        Specifies the output path. Resulting trie will be
                        written as *.trie file.
  -b BLOSUM, --blosum BLOSUM
                        Specifies BLOSUM distance matrix (default BLOSUM50;
                        available BLOSUM45, BLOSUM90).

"""

import sys
import argparse
import csv
import logging
import os
import pandas
import itertools as itr

from Distance2SelfBinding import Distance2Self
from DistanceMatrix import DistanceMatrix

def load_blossum(blos):
    """
    loads a BLOSUM matrix

    :param str blos: Specifeis the BLOSUM matrix to lead
    :return: dict(str, dict(str, float)) - A BLOSUM1 matrix
    """
    try:
        mod = __import__('DistanceMatrices', fromlist=[blos])
        return DistanceMatrix(getattr(mod, blos))
    except:
        mod = __import__('DistanceMatrices', fromlist=["BLOSUM50_distances"])
        return DistanceMatrix(getattr(mod, "BLOSUM50_distances"))

def main():
    parser = argparse.ArgumentParser(
        description="Generation of trie data structure for a set of proteins for distance to self calculation",

    )

    parser.add_argument("-i", "--input",
                        required=True,
                        type=str,
                        help="Protein sequences in fasta file format. All peptides of the specified length will be generated for the trie generation.",
                        )

    parser.add_argument("-l", "--length",
                             required=False,
                             default=9,
                             type=int,
                             help="Specifies the length of peptides which will be generated."
                             )

    parser.add_argument("-o", "--output",
                        required=True,
                        type=str,
                        help="Specifies the output path. Resulting trie will be written as *.trie file.",
                        )

    parser.add_argument("-b", "--blosum",
                        required=False,
                        default="BLOSUM50",
                        type=str,
                        help="Specifies BLOSUM distance matrix (default BLOSUM50; available BLOSUM45, BLOSUM90).", 
                        )

    args = parser.parse_args()
    blos = load_blossum("{blos}_distances".format(blos=args.blosum.strip().upper()))
    dist2self = Distance2Self(blos,saveTrieFile=True)

    dist2self.generate_trie(args.input, peptideLength=args.length, outfile=args.output, peplength=args.length)

if __name__ == "__main__":
    sys.exit(main())