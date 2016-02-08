#!/usr/bin/env python
"""
Commandline tool for tap prediction

usage: taprediction.py [-h] [-m {svmtap,smmtap,doytchinova}] [-v VERSION] -i
                       INPUT [-t {fasta,peptide}] [-l LENGTH] [-op OPTIONS] -o
                       OUTPUT

Commandline tool for TAP prediction

optional arguments:
  -h, --help            show this help message and exit
  -m {svmtap,smmtap,doytchinova}, --method {svmtap,smmtap,doytchinova}
                        The name of the prediction method
  -v VERSION, --version VERSION
                        The version of the prediction method
  -i INPUT, --input INPUT
                        Path to the input file
  -t {fasta,peptide}, --type {fasta,peptide}
                        The data type of the input (fasta, peptide list)
  -l LENGTH, --length LENGTH
                        The length of peptides
  -op OPTIONS, --options OPTIONS
                        Additional options that get directly past to the tool
  -o OUTPUT, --output OUTPUT
                        Path to the output file

"""
import sys
import argparse

from Fred2.Core import Protein, Peptide, Allele
from Fred2.IO import read_fasta
from Fred2.TAPPrediction import TAPPredictorFactory
from Fred2.Core import generate_peptides_from_proteins

def read_lines(file, in_type=Peptide):
    peptides = []
    with open(file, "r") as f:
        for l in f:
            if not l.startswith("#") and l.strip() != "" and not l.startswith("Epitope") and not l.startswith(
                    "Sequence"):
                print l, l.split()
                pep = l.split()[0].strip()
                peptides.append(in_type(pep))
    return peptides

def main():
    model = argparse.ArgumentParser(
        description='Commandline tool for TAP prediction',
        )

    model.add_argument('-m',
        '--method',
        type=str,
        choices=TAPPredictorFactory.available_methods().keys(),
        default="svmtap",
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

    model.add_argument('-t',
        '--type',
        choices=["fasta", "peptide"],
        type=str,
        default="fasta",
        help='The data type of the input (fasta, peptide list)'
        )

    model.add_argument('-l',
        '--length',
        type=int,
        default=9,
        help='The length of peptides'
        )

    model.add_argument('-op',
        '--options',
        type=str,
        default="",
        help="Additional options that get directly past to the tool"
    )

    model.add_argument('-o',
        '--output',
        type=str,
        required=True,
        help='Path to the output file'
        )

    args = model.parse_args()

    #fasta protein
    if args.type == "fasta":
        with open(args.input, 'r') as f:
            first_line = f.readline()
        sep_pos = 1 if first_line.count("|") else 0
        proteins = read_fasta(args.input, in_type=Protein, id_position=sep_pos)
        peptides = generate_peptides_from_proteins(proteins, int(args.length))
    elif args.type == "peptide":
        peptides = read_lines(args.input, in_type=Peptide)
    else:
        sys.stderr.write('Input type not known\n')
        return -1

    if args.version == "":
        result = TAPPredictorFactory(args.method).predict(peptides, options=args.options)
    else:
        result = TAPPredictorFactory(args.method, version=args.version).predict(peptides, options=args.options)

    #write to TSV columns sequence method score...,protein-id/transcript-id
    with open(args.output, "w") as f:
        proteins = "\tProtein ID" if args.type == "fasta" else ""
        f.write("Sequence\tMethod\t"+"Score"+proteins+"\n")
        for index, row in result.iterrows():
            p = index
            proteins = ",".join(prot.transcript_id for prot in p.get_all_proteins()) if args.type == "fasta" else ""
            f.write(str(p)+"\t"+"\t".join("%s\t%.3f"%(method, score) for
                                          method, score in row.iteritems())+"\t"+proteins+"\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())