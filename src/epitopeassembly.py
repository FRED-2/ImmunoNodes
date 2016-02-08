#!/usr/bin/env python
"""
Epitope assembly command line tool

usage: epitopeassembly.py [-h] -i INPUT
                          [-m {proteasmm_i,netchop,proteasmm_c,pcm}]
                          [-v VERSION] [-w WEIGHT] [-a] [-s SOLVER]
                          [-so SOLVER_OPTIONS] -o OUTPUT

Commandline tool for string-of-beads epitope assembly

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the input file
  -m {proteasmm_i,netchop,proteasmm_c,pcm}, --method {proteasmm_i,netchop,proteasmm_c,pcm}
                        The name of the prediction method
  -v VERSION, --version VERSION
                        The version of the prediction method
  -w WEIGHT, --weight WEIGHT
                        Specifies how strong unwanted cleavage sites should be
                        punished [0,1], where 0 means they will be ignored,
                        and 1 the sum of all unwanted cleave sites is
                        subtracted from the cleave site between two epitopes
  -a, --approximate     Flag whether epitope ordering should be approximated
  -s SOLVER, --solver SOLVER
                        The ILP solver to be used by Pyomo (must be in PATH)
  -so SOLVER_OPTIONS, --solver_options SOLVER_OPTIONS
                        Solver specific options (will not be checked for
                        correctness)
  -o OUTPUT, --output OUTPUT
                        Path to the output file

"""
import sys
import argparse

from Fred2.Core import Peptide
from Fred2.CleavagePrediction import CleavageSitePredictorFactory
from Fred2.EpitopeAssembly import EpitopeAssembly

def read_lines(file, in_type=Peptide):
    peptides = []

    with open(file, "r") as f:
        for l in f:
            if not l.startswith("#") and l.strip() != "" and not l.startswith("Epitope") and not l.startswith("Sequence"):
                pep = l.split()[0].strip()
                peptides.append(in_type(pep))
    return peptides

def main():
    model = argparse.ArgumentParser(
        description='Commandline tool for string-of-beads epitope assembly',
        )

    model.add_argument('-i',
        '--input',
        type=str,
        required=True,
        help='Path to the input file'
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

    model.add_argument('-w',
        '--weight',
        type=float,
        default=0.0,
        help="Specifies how strong unwanted cleavage sites should be punished [0,1], \
                             where 0 means they will be ignored, and 1 the sum of all unwanted cleave sites is \
                             subtracted from the cleave site between two epitopes"
        )

    model.add_argument('-a',
        '--approximate',
        action="store_true",
        help="Flag whether epitope ordering should be approximated"
    )

    model.add_argument('-s',
        '--solver',
        type=str,
        default="cbc",
        help="The ILP solver to be used by Pyomo (must be in PATH)"
    )

    model.add_argument('-so',
        '--solver_options',
        type=str,
        default="",
        help="Solver specific options (will not be checked for correctness)"
    )

    model.add_argument('-o',
        '--output',
        type=str,
        required=True,
        help='Path to the output file'
        )

    args = model.parse_args()

    try:
        peptides = read_lines(args.input)
    except Exception as e:
        print e
        sys.stderr.write("Input file could not be read. Please check the existence and input format.")
        return -1

    try:
        if args.version == "":
            predictor = CleavageSitePredictorFactory(args.method)
        else:
            predictor = CleavageSitePredictorFactory(args.method, version=args.version)
    except ValueError as e:
        sys.stderr.write(str(e))
        return -1

    try:
        assembler = EpitopeAssembly(peptides, predictor, solver=args.solver, weight=args.weight)
    except Exception as e:
        sys.stderr.write(str(e))
        return -1

    try:
        if args.approximate:
            assembly = assembler.approximate()
        else:
            options = {} 
            for opt in args.solver_options.split():
                name,value = opt.partition("=")[::2]
                try:
                    options[name] = float(value)
                except Exception:
                    options[name] = value
            assembly = assembler.solve(options=options)
    except Exception as e:
        sys.stderr.write(str(e))
        return -1

    try:
        with open(args.output, "w") as f:
            f.write(">assembled_polypeptide\n"+"".join(str(p) for p in assembly))
    except IOError as e:
        sys.stderr.write(str(e))
        return -1

    return 0

if __name__ == "__main__":
    sys.exit(main())