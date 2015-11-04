#!/usr/bin/env python
"""
Epitope assembly command line tool

"""
import sys

from CTDopts.CTDopts import CTDModel

from Fred2.Core import Peptide
from Fred2.IO import read_lines
from Fred2.CleavagePrediction import CleavageSitePredictorFactory
from Fred2.EpitopeAssembly import EpitopeAssembly



def main():
#Specify CTD interface
    # Every CTD Model has to have at least a name and a version, plus any of the optional attributes below them.
    model = CTDModel(
        name='EpitopeAssembly',  # required
        version='1.0',  # required
        description='Commandline tool for string-of-beads epitope assembly',
        manual='manual string',
        executableName='epitopeassembly',
        )

    model.add(
        'input',
        type="input-file",
        required=True,
        description='Path to the input file'
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
        'weight',
        type=float,
        bounds=(0, None),
        default=0,
        description="Specifies how strong unwanted cleavage sites should be punished [0,1], \
                             where 0 means they will be ignored, and 1 the sum of all unwanted cleave sites is \
                             subtracted from the cleave site between two epitopes"
        )

    model.add(
        'approximate',
        type=bool,
        default=False,
        description="Flag whether epitope ordering should be approximated"
    )

    model.add(
        'solver',
        type=str,
        default="cbc",
        description="The ILP solver to be used by Pyomo (must be in PATH)"
    )

    model.add(
        'solver-options',
        type=str,
        default=None,
        description="Solver specific options (will not be checked for correctness)"
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

    try:
        peptides = read_lines(args["input"], type=Peptide)
    except Exception as e:
        sys.stderr.write("Input file could not be read. Please check the existence and input format.")
        return -1

    try:
        if args["version"] is None:
            predictor = CleavageSitePredictorFactory(args["method"])
        else:
            predictor = CleavageSitePredictorFactory(args["method"], version=args["version"])
    except ValueError as e:
        sys.stderr.write(str(e))
        return -1

    try:
        assembler = EpitopeAssembly(peptides, predictor, solver=args["solver"], weight=args["weight"])
    except Exception as e:
        sys.stderr.write(str(e))
        return -1

    try:
        if args["approximate"]:
            assembly = assembler.approximate()
        else:
            assembly = assembler.solve(options=args["solver-options"])
    except Exception as e:
        sys.stderr.write(str(e))
        return -1

    try:
        with open(args["output"], "w") as f:
            f.write(">assembled_polypeptide\n"+"".join(str(p) for p in assembly))
    except IOError as e:
        sys.stderr.write(str(e))
        return -1

    return 0

if __name__ == "__main__":
    sys.exit(main())