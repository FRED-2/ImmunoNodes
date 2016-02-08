#!/usr/bin/env python
"""
Cleavage prediction src tool

usage: cleavageprediction.py [-h] [-m {proteasmm_i,netchop,proteasmm_c,pcm}]
                             [-v VERSION] -i INPUT [-l LENGTH] [-op OPTIONS]
                             -o OUTPUT

Commandline tool for cleavage site prediction

optional arguments:
  -h, --help            show this help message and exit
  -m, --method {proteasmm_i,netchop,proteasmm_c,pcm}
                        The name of the prediction method
  -v VERSION, --version VERSION
                        The version of the prediction method
  -i INPUT, --input INPUT
                        Path to the input file (in fasta format)
  -l LENGTH, --length LENGTH
                        The length of peptides
  -op OPTIONS, --options OPTIONS
                        Additional options that get directly past to the tool
  -o OUTPUT, --output OUTPUT
                        Path to the output file

"""
import sys
import argparse

from Fred2.Core import Protein
from Fred2.IO import read_fasta
from Fred2.CleavagePrediction import CleavageSitePredictorFactory

def main():
    #Specify CTD interface
    # Every CTD Model has to have at least a name and a version, plus any of the optional attributes below them.
    model = argparse.ArgumentParser(
        description='Commandline tool for cleavage site prediction',
        )

    model.add_argument('-m',
        '--method',
        type=str,
        choices=CleavageSitePredictorFactory.available_methods().keys(),
        default="pcm",
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
        help='Path to the input file (in fasta format)'
        )

    model.add_argument('-l',
        '--length',
        type=int,
        default=0,
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
    peptides = read_fasta(args.input, in_type=Protein)

    if args.version == "":
        predictor = CleavageSitePredictorFactory(args.method)
        result = predictor.predict(peptides, options=args.method)
    else:
        predictor = CleavageSitePredictorFactory(args.method, version=args.version)
        result = predictor.predict(peptides, options=args.method)

    #if length is specified, than generate compact output
    if int(args.length) > 0:
        length = int(args.length)
        with open(args.output, "w") as f:
            f.write("Sequence\tMethod\tScore\tProtein ID\tPosition\n")
            for seq_id in set(result.index.get_level_values(0)):
                    seq = "".join(result.ix[seq_id]["Seq"])
                    for start in xrange(len(seq)-(length-1)):
                        pep_seq = seq[start:(start+length)]
                        score = result.loc[(seq_id, start+(length-1)), predictor.name]
                        f.write(pep_seq+"\t"+predictor.name+"\t"+"%.3f"%score+"\t"+seq_id+"\t"+str(start)+"\n")
    else:
        result.to_csv(args.output, float_format="%.3f", sep="\t")
    return 0

if __name__ == "__main__":
    sys.exit(main())