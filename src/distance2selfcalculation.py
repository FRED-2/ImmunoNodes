#!/usr/bin/env python
"""
Command line tool for distance to self calculation

usage: distance2selfcalculation.py [-h] -i INPUT [-s SEQUENCE]
                                   [--precalculated_trie {uniprot_proteome_l8, uniprot_proteome_l9, uniprot_proteome_l10, uniprot_proteome_l11} | --custom_trie CUSTOM_TRIE]
                                   [-k K] [-b BLOSUM] -o OUTPUT

Distance to self calculation

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Tab-separated peptide file (e.g. from
                        epitopeprediction)
  -s SEQUENCE, --sequence SEQUENCE
                        The columns name of the peptide sequences
  --precalculated_trie {uniprot_proteome_l8,uniprot_proteome_l9,uniprot_proteome_l10,uniprot_proteome_l11}
                        Use one of the pre-generated tries (default
                        uniprot_proteome_l9; available uniprot_proteome_l8,uniprot_proteome_l10,
                        uniprot_proteome_l11). Pre-generated tries have been generated
                        from uniprot human proteome for different peptide
                        lengths. May not be combined with parameter
                        --custom_trie
  --custom_trie CUSTOM_TRIE
                        Path to a self-generated custom trie. May not be
                        combined with parameter --precalculated_trie
  -k K, --k K           Specifies the number of closest self-peptides to find
  -b BLOSUM, --blosum BLOSUM
                        Specifies BLOSUM distance matrix (default BLOSUM50;
                        available BLOSUM45, BLOSUM90)
  -o OUTPUT, --output OUTPUT
                        Specifies the output path. Results will be written to
                        CSV

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
    :return: dict(str, dict(str, float)) - A BLOSUM matrix
    """
    try:
        mod = __import__('DistanceMatrices', fromlist=[blos])
        return DistanceMatrix(getattr(mod, blos))
    except:
        mod = __import__('DistanceMatrices', fromlist=["BLOSUM50_distances"])
        return DistanceMatrix(getattr(mod, "BLOSUM50_distances"))

def main():
    parser = argparse.ArgumentParser(
        description="Distance to self calculation",

    )

    parser.add_argument("-i", "--input",
                        required=True,
                        type=str,
                        help="Tab-separated peptide file (e.g. from epitopeprediction)",
                        )

    parser.add_argument("-s", "--sequence",
                        required=False,
                        default="neopeptide",
                        type=str,
                        help="The columns name of the peptide sequences",
                        )

    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument("--precalculated_trie",
                        default="uniprot_proteome_l9",
                        choices=["uniprot_proteome_l8","uniprot_proteome_l9","uniprot_proteome_l10","uniprot_proteome_l11"],
                        type=str,
                        help="Use one of the pre-generated tries (default uniprot_proteome_l9; available uniprot_proteome_l8, uniprot_proteome_l10, uniprot_proteome_l11). Pre-generated tries have been generated from uniprot human proteome for different peptide lengths. May not be combined with parameter --custom_trie")
    group.add_argument('--custom_trie',
                        type=str,
                        help="Path to a self-generated custom trie. May not be combined with parameter --precalculated_trie")
    parser.add_argument("-k", "--k",
                             required=False,
                             default=1,
                             type=int,
                             help="Specifies the number of closest self-peptides to find"
                             )
    parser.add_argument("-b", "--blosum",
                        required=False,
                        default="BLOSUM50",
                        type=str,
                        help="Specifies BLOSUM distance matrix (default BLOSUM50; available BLOSUM45, BLOSUM90)",
                        )
    parser.add_argument("-o", "--output",
                        required=True,
                        type=str,
                        help="Specifies the output path. Results will be written to CSV",
                        )

    args = parser.parse_args()
    blos = load_blossum("{blos}_distances".format(blos=args.blosum.strip().upper()))
    dist2self = Distance2Self(blos,saveTrieFile=True)
    df = pandas.read_csv(args.input, sep="\t", index_col=False)
    peps = list(set(df[args.sequence]))

    if args.custom_trie is not None:
        specifiedTrie = args.custom_trie
    else:
        specifiedTrie = args.precalculated_trie

    res = dist2self.calculate_distances(peps, pep_header=args.sequence, specifiedTrie=specifiedTrie, n=args.k)
    merged = pandas.merge(df, res, how="outer",on=[args.sequence])
    merged.to_csv(args.output, sep="\t",index=False)

if __name__ == "__main__":
    sys.exit(main())