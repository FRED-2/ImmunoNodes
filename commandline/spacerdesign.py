#!/usr/bin/env python
# coding=utf-8
"""


###################################################################
Designing String-of-beads with optimal spacers
###################################################################
Authors: Benjamin Schubert and Oliver Kohlbacher
Date: June 2015
Version: 1.0
License: This software is under a three-clause BSD license
Introduction:
-------------
The software is a novel approach to construct epitope-based string-of-beads
vaccines in optimal order and with sequence-optimized spacers of flexible length
such that the recovery of contained epitopes is maximized and immunogenicity of 
arising neo-epitopes is reduced. 
Requirement:
-------------
Spacer Design uses the following software and libraries:
1) Python 2.7 (https://www.python.org/)
2) Fred2 (https://github.com/FRED-2/Fred2)
3) Cplex >= 12.5 (www.ilog.com)
4) LKH TSP-Approximation >= 2.0.7 (http://www.akira.ruc.dk/~keld/research/LKH/)
Please make sure you have installed said software/libraries
and their dependencies.
Installation:
-------------
First install all required software and libraries. CPLEX/LKH should be globally executable
via command line. 
Usage:
-------------
usage: OptiVac.py [-h] -i INPUT -a ALLELES [-k MAX_LENGTH] [-al ALPHA]
                   [-be BETA] [-cp CLEAVAGE_PREDICTION]
                   [-ep EPITOPE_PREDICTION] [-thr THRESHOLD] -o OUTPUT
                   [-t THREADS]
The software is a novel approach to construct epitope-based string-of-beads
vaccines in optimal order and with sequence-optimized spacers of flexible
length such that the recovery of contained epitopes is maximized and
immunogenicity of arising neo-epitopes is reduced.
optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        File containing epitopes (one peptide per line)
  -a ALLELES, --alleles ALLELES
                        Specifies file containing HLA alleles with
                        corresponding HLA probabilities (one HLA per line)
  -k MAX_LENGTH, --max_length MAX_LENGTH
                        Specifies the max. length of the spacers (default 6)
  -al ALPHA, --alpha ALPHA
                        Specifies the first-order preference of the user in
                        the model [0,1] (default 0.99)
  -be BETA, --beta BETA
                        Specifies the second-order preference of the user in
                        the model [0,1] (default 0).
  -cp CLEAVAGE_PREDICTION, --cleavage_prediction CLEAVAGE_PREDICTION
                        Specifies the used cleavage prediction method (default
                        PCM) [available: PCM, ProteaSMMConsecutive, ProteaSMMImmuno]
  -ep EPITOPE_PREDICTION, --epitope_prediction EPITOPE_PREDICTION
                        Specifies the used epitope prediction method (default
                        Syfpeithi) [available: Syfpeithi, BIMAS, SMM, SMMPMBEC]
  -thr THRESHOLD, --threshold THRESHOLD
                        Specifies epitope prediction threshold for SYFPEITHI
                        (default 20).
  -o OUTPUT, --output OUTPUT
                        Specifies the output file.
  -t THREADS, --threads THREADS
                        Specifies number of threads. If not specified all
                        available logical cpus are used.
"""

import sys
import multiprocessing as mp

from Fred2.IO import FileReader
from Fred2.Core import Allele
from Fred2.Core import Peptide
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.EpitopeAssembly.EpitopeAssembly import EpitopeAssemblyWithSpacer
from Fred2.CleavagePrediction import CleavageSitePredictorFactory

from CTDopts.CTDopts import CTDModel


def generate_alleles(allele_file, generated=None):
    """
                generate allele objects from input
    """
    result=[]
    with open(allele_file, "r") as f:
        for l in f:
            al, freq = l.replace(",", " ").replace(";", " ").replace("\n", "").split()
            if al.split("HLA-")[-1][0] in ["A","B","C"]:
                result.append(Allele(al, prob=float(freq)))
    return result


def main():
    parser = CTDModel(
        name='SpacerDesign',  # required
        version='1.0',  # required
        description='he software is a novel approach to construct epitope-based string-of-beads \
vaccines in optimal order and with sequence-optimized spacers of flexible length \
such that the recovery of contained epitopes is maximized and immunogenicity of \
arising neo-epitopes is reduced.',
        manual='manual string',
        executableName='spacerdesign',
        )

    parser.add("input",
                        required=True,
                        help="File containing epitopes (one peptide per line)",
                        short_name="i"
    )

    parser.add("alleles",
                        required=True,
                        help="Specifies file containing HLA alleles with corresponding HLA probabilities (one HLA per line)",
                        short_name="a"

    )

    #parameters of the model
    parser.add("max_length",
                        default=6,
                        type=int,
                        help="Specifies the max. length of the spacers (default 6)",
                        short_name="k")

    parser.add("alpha",
                        default=0.99,
                        type=float,
                        num_range=(0.0, 1.0),
                        help="Specifies the first-order preference of the user in the model [0,1] (default 0.99)",
                        short_name="al")
    parser.add("beta",
                        default=0.0,
                        type=float,
                        num_range=(0.0, 1.0),
                        help="Specifies the second-order preference of the user in the model [0,1] (default 0).",
                        short_name="be")

    parser.add("cleavage-prediction",
                        default="PCM",
                        choice=["PCM", "PROTEASMM_C", "PROTEASMM_i"],
                        help="Specifies the used cleavage prediction method (default PCM) [available: PCM, PROTEASMM_C, PROTEASMM_S]",
                        short_name="cp"
    )
    parser.add("cleavage-version",
                        default=None,
                        help="Specifies the ",
                        short_name="cp"
    )
    parser.add("epitope-prediction",
                        default="Syfpeithi",
                        choice=["Syfpeithi", "SMM", "SMMPMBEC", "BIMAS"],
                        help="Specifies the used epitope prediction method (default Syfpeithi) [available: Syfpeithi, BIMAS, SMM, SMMPMBEC]",
                        short_name="ep"
    )
    parser.add("threshold",
                        default=20,
                        type=float,
                        help="Specifies epitope prediction threshold for SYFPEITHI (default 20).",
                        short_name="thr")

    parser.add("output",
                        required=True,
                        help="Specifies the output file.",
                        short_name="o")

    parser.add("threads",
                        type=int,
                        default=None,
                        help="Specifies number of threads. If not specified all available logical cpus are used.",
                        short_name="t")


    parser.add(
        'ctdout',
        default=None,
        type="output-file",
        description='Output path to for cds'
        )

    args_str = sys.argv[1:] if sys.argv[1:] else ["--help"]
    args = parser.parse_cl_args(cl_args=args_str)

    if args["ctdout"] is not None:
        parser.write_ctd(args[args["ctdout"]])
        return 0

    #parse input
    peptides = list(FileReader.read_lines(args["input"], type=Peptide))
    #read in alleles
    alleles = generate_alleles(args["alleles"])

    if args["cleavage-prediction"].upper() not in ["PCM", "PROTEASMM_C", "PROTEASMM_S"]:
        sys.stderr.write("Specified cleavage predictor is currently not supported. \
                         Please choose either PCM, PROTEASMM_C, or PROTEASMM_S")
        sys.exit(-1)

    if args["epitope-prediction"].upper() not in ["SYFPEITHI", "BIMAS", "SMM", "SMMPMBEC"]:
        sys.stderr.write("Specified cleavage predictor is currently not supported. \
                         Please choose either Syfpeithi, BIMAS, SMM, SMMPMBEC")
        sys.exit(-1)

    #set-up model
    cl_pred = CleavageSitePredictorFactory(args["cleavage-prediction"], version=args["cleavage-version"])
    epi_pred = EpitopePredictorFactory(args["epitope_prediction"], version=args["epitope-version"])


    thr = {a.name:args.threshold for a in alleles}

    solver = EpitopeAssemblyWithSpacer(peptides,cl_pred,epi_pred,alleles,
                                       k=args.max_length,en=9,threshold=thr,
                                       solver="cbc", alpha=args.alpha, beta=args.beta,
                                       verbosity=0)

    #solve
    #pre-processing has to be disable otherwise many solver will destroy the symmetry of the problem
    #how to do this is dependent on the solver used. For CPLEX it is preprocessing_presolve=n
    #TODO:CBC should be shipped with the node
    #TODO: has to be tested with CBC
    #TODO: LHK has to be shipped as well -> only academic license!
    threads = mp.cpu_count() if args.threads is None else args.threads
    svbws = solver.approximate(threads=threads,options="preprocess=off threads=1")

    with open(args.output, "w") as f:
        f.write(">assembled_spacer_design")
        f.write("".join(map(str,svbws)))
    return 0

if __name__ == "__main__":
    sys.exit(main())