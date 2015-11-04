#!/usr/bin/env python
"""
Commandline tool for epitope prediction


"""
import sys

from CTDopts.CTDopts import CTDModel

from Fred2.Core import Protein, Peptide, Allele
from Fred2.IO import read_lines, read_fasta
from Fred2.EpitopePrediction import EpitopePredictorFactory
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
        choices=EpitopePredictorFactory.available_methods().keys(),
        default="bimas",
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
        type="input-file",
        required=True,
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
        num_range=(8, 16),
        type=int,
        default=9,
        description='The length of peptides'
        )

    model.add(
        'alleles',
        type="input-file",
        required=True,
        description='Path to the allele file (one per line in new nomenclature)'
        )

    model.add(
        'options',
        type=str,
        default=None,
        description="Additional options that get directly past to the tool"
    )

    model.add(
        'output',
        type="output-file",
        required=True,
        description='Path to the output file'
        )

    model.add(
        'ctdout',
        default=None,
        type="output-file",
        description='Output path to for cds'
        )

    args_str = sys.argv[1:] if sys.argv[1:] else ["--help"]
    args = model.parse_cl_args(cl_args=args_str)

    if args["ctdout"] is not None:
        model.write_ctd(args[args["ctdout"]])
        return 0

    #fasta protein
    if args["type"] == "fasta":
        proteins = read_fasta(args["input"], type=Protein)
        peptides = generate_peptides_from_protein(proteins, args["length"])
    elif args["type"] == "peptide":
        peptides = read_lines(args["input"], type=Peptide)
    else:
        sys.stderr.write('Input type not known\n')
        return -1

    #read in alleles
    alleles = read_lines(args["alleles"], type=Allele)
    if args["version"] is None:
        result = EpitopePredictorFactory(args["method"]).predict(peptides, alleles, options=args["options"])
    else:
        result = EpitopePredictorFactory(args["method"], version=args["version"]).predict(peptides, alleles,
                                                                 options=args["options"])
    print result
    result.to_csv(args["output"])
    return 0

if __name__ == "__main__":
    sys.exit(main())