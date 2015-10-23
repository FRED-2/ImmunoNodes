#!/usr/bin/env python
"""
Cleavage prediction commandline tool

"""
import sys
import pandas

from CTDopts.CTDopts import CTDModel

from Fred2.Core import Protein, Peptide, Allele
from Fred2.IO import read_lines, read_fasta
from Fred2.CleavagePrediction import CleavageSitePredictorFactory, CleavageFragmentPredictionResult
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
        choices=CleavageSitePredictorFactory.available_methods().keys(),
        default="pcm",
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
        num_range=(1, None),
        type=int,
        default=None,
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
        type="output-file",
        description='Path to the output file'
        )
    args_str = sys.argv[1:] if sys.argv[1:] else ["--help"]
    args = model.parse_cl_args(cl_args=args_str)

    #fasta protein
    if args["type"] == "fasta":
        peptides = read_fasta(args["input"], type=Protein)
        if args["length"] is not None:
            peptides = generate_peptides_from_protein(peptides, args["length"])
    elif args["type"] == "peptide":
        peptides = read_lines(args["input"], type=Peptide)
    else:
        sys.stderr.write('Input type not known\n')
        return -1

    if args["version"] is None:
        predictor = CleavageSitePredictorFactory(args["method"]).predict(peptides, options=args["options"])
        result = predictor.predict(peptides, options=args["options"])
    else:
        predictor = CleavageSitePredictorFactory(args["method"], version=args["version"])
        result = predictor.predict(peptides, options=args["options"])

    #if length is specified, than generate compact output
    if args["length"] is not None:
        length = args["length"]
        r_comp = {predictor.name:{}}
        for seq_id in set(result.index.get_level_values(0)):
                seq = "".join(result.ix[seq_id]["Seq"])
                for start in xrange(len(seq)-(length-1)):
                    pep_seq = seq[start:(start+(length-1))]
                    r_comp[predictor.name][pep_seq] = result.loc[(seq_id++(length-1), start), predictor.name]
        result = CleavageFragmentPredictionResult.from_dict(result)
        result.index = pandas.MultiIndex.from_tuples([tuple((i, predictor.name)) for i in result.index],
                                                        names=['Seq', 'Method'])
    result.to_csv(args["out"])
    return 0