#!/usr/bin/env python
"""
Command line tool to assign allele frequencies to alleles lists

usage: allelefrequency.py [-h] -p
                          {Oceania,South-East_Asia,South-West_Asia,Europe,North_Africa,Australia,South_America,Sub-Saharan_Africa,Other,North-East_Asia,North_America}
                          [-a ALLELES] [-t THRESHOLD] -o OUTPUT

Commandline tool for allele frequency assignment

optional arguments:
  -h, --help            show this help message and exit
  -p {Oceania,South-East_Asia,South-West_Asia,Europe,North_Africa,Australia,South_America,Sub-Saharan_Africa,Other,North-East_Asia,North_America}, --population {Oceania,South-East_Asia,South-West_Asia,Europe,North_Africa,Australia,South_America,Sub-Saharan_Africa,Other,North-East_Asia,North_America}
                        Specifies the population (e.g. Europe)
  -a ALLELES, --alleles ALLELES
                        Path to the allele file (one per line in new
                        nomenclature)
  -t THRESHOLD, --threshold THRESHOLD
                        Frequency threshold to include alleles
  -o OUTPUT, --output OUTPUT
                        Path to the output file

"""
from Fred2.Core import Allele
from Fred2.IO import read_lines

import sys
from data.geo import geo
import argparse


def main():
    model = argparse.ArgumentParser(
        description='Commandline tool for allele frequency assignment',
        )

    model.add_argument(
        '-p','--population',
        choices=geo.keys(),
        type=str,
        required=True,
        help='Specifies the population (e.g. Europe)',
        )

    model.add_argument(
        '-a','--alleles',
        type=str,
        default="",
        help='Path to the allele file (one per line in new nomenclature)',
        )

    model.add_argument(
        '-t','--threshold',
        type=float,
        default=0.001,
        help="Frequency threshold to include alleles",
    )

    model.add_argument(
        '-o','--output',
        type=str,
        required=True,
        help='Path to the output file',
        )

    args = model.parse_args()

    thr = float(args.threshold)
    with open(args.output, "w") as f:
        if args.population in geo:
            freqs = geo[args.population]
        else:
            sys.stderr("{pop} could not be found".format(pop=args.population))
            return -1

        if args.alleles != "":
            alleles = read_lines(args.alleles, in_type=Allele)
            for a in alleles:
                fr = freqs.get(a.name, 0.0)
                if fr >= thr:
                    f.write(a.name+"\t"+"%.3f"%fr+"\n")
        else:
            for a, fr in freqs.iteritems():
                fr = fr
                if fr >= thr:
                    f.write(a+"\t"+"%.3f"%fr+"\n")
    return 0

if __name__ == "__main__":
    sys.exit(main())