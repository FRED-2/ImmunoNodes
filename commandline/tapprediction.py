#!/usr/bin/env python
"""
Commandline tool for tap prediction


"""
import sys

from CTDopts import CTDModel

from Fred2.Core import Protein, Peptide, Allele
from Fred2.IO import read_lines, read_fasta
from Fred2.TAPPrediction import TAPPredictorFactory
from Fred2.Core import generate_peptides_from_protein


def main():
#Specify CTD interface
    # Every CTD Model has to have at least a name and a version, plus any of the optional attributes below them.
    model = CTDModel(
        name='EpitopePredicton',  # required
        version='1.0',  # required
        description='Commandline tool for epitope prediction',
        manual='manual string',
        executableName='epitopeprediction',
        )

    model.add(
        'method',
        type=str,
        choices=TAPPredictorFactory.available_methods().keys(),
        default="svmtap",
        description='The name of the prediction method'
        )

    model.add(
        'version',
        type=str,
        default=None,
        description='The version of the prediction method'
        )

    model.add(
        'input',
        type="input_file",
        description='Path to the input file'
        )

    model.add(
        'type',
        choices=["fasta","peptide"],
        type=str,
        default="fasta",
        description='The data type of the input (fasta, peptide list)'
        )

    model.add(
        'length',
        num_range=(9, 16),
        type=int,
        default=9,
        description='The length of peptides'
        )

    model.add(
        'options',
        type=str,
        default=None,
        description="Additional options that get directly past to the tool"
    )

    model.add(
        'output',
        type="output_file",
        description='Path to the output file'
        )

    args = model.parse_cl_args(cl_args=sys.argv[1:])

    #fasta protein
    if args["type"] == "fasta":
        proteins = read_fasta(args["input"], type=Protein)
        peptides = generate_peptides_from_protein(proteins, args["length"])
    elif args["type"] == "peptide":
        peptides = read_lines(args["input"], type=Peptide)
    else:
        sys.stderr.write('Input type not known\n')
        return -1

    if args["version"] is None:
        result = TAPPredictorFactory(args["method"]).predict(peptides, options=args["options"])
    else:
        result = TAPPredictorFactory(args["method"], version=args["version"]).predict(peptides, options=args["options"])
    result.to_csv(args["out"])
    return 0


if __name__ == "__main__":
    sys.exit(main())