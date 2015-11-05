#!/usr/bin/env python
"""
Cleavage prediction commandline tool

"""
import sys

from CTDopts.CTDopts import CTDModel

from Fred2.Core import Protein
from Fred2.IO import read_fasta
from Fred2.CleavagePrediction import CleavageSitePredictorFactory


def main():
#Specify CTD interface
    # Every CTD Model has to have at least a name and a version, plus any of the optional attributes below them.
    model = CTDModel(
        name='CleavagePrediction',  # required
        version='1.0',  # required
        description='Commandline tool for cleavage site prediction',
        manual='manual string',
        executableName='cleavageprediction',
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
        required=True,
        description='Path to the input file (in fasta format)'
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
    peptides = read_fasta(args["input"], type=Protein)

    if args["version"] is None:
        predictor = CleavageSitePredictorFactory(args["method"])
        result = predictor.predict(peptides, options=args["options"])
    else:
        predictor = CleavageSitePredictorFactory(args["method"], version=args["version"])
        result = predictor.predict(peptides, options=args["options"])

    #if length is specified, than generate compact output
    if args["length"] is not None:
        length = int(args["length"])
        with open(args["output"], "w") as f:
            f.write("Sequence\tMethod\tScore\tProtein ID\tPosition\n")
            for seq_id in set(result.index.get_level_values(0)):
                    seq = "".join(result.ix[seq_id]["Seq"])
                    for start in xrange(len(seq)-(length-1)):
                        pep_seq = seq[start:(start+(length-1))]
                        score = result.loc[(seq_id, start+(length-2)), predictor.name]
                        f.write(pep_seq+"\t"+predictor.name+"\t"+"%.3f"%score+"\t"+seq_id+"\t"+str(start)+"\n")
    else:
        result.to_csv(args["output"], float_format="%.3f", sep="\t")
    return 0

if __name__ == "__main__":
    sys.exit(main())