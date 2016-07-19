#!/usr/bin/env python
"""
Commandline tool for epitope prediction

usage: epitopeprediction.py [-h]
                            [-m {netmhc,smmpmbec,syfpeithi,netmhcpan,netctlpan,smm,tepitopepan,arb,pickpocket,epidemix,netmhcii,netmhciipan,comblibsidney,unitope,hammer,svmhc,bimas}]
                            [-v VERSION] -i INPUT [-t {fasta,peptide}]
                            [-l {8,9,10,11,12,13,14,15,16,17}] -a ALLELES
                            [-op OPTIONS] -o OUTPUT

Process some integers.

optional arguments:
  -h, --help            show this help message and exit
  -m {netmhc,smmpmbec,syfpeithi,netmhcpan,netctlpan,smm,tepitopepan,arb,pickpocket,epidemix,netmhcii,netmhciipan,comblibsidney,unitope,hammer,svmhc,bimas}, --method {netmhc,smmpmbec,syfpeithi,netmhcpan,netctlpan,smm,tepitopepan,arb,pickpocket,epidemix,netmhcii,netmhciipan,comblibsidney,unitope,hammer,svmhc,bimas}
                        The name of the prediction method
  -v VERSION, --version VERSION
                        The version of the prediction method
  -i INPUT, --input INPUT
                        Path to the input file
  -t {fasta,peptide}, --type {fasta,peptide}
                        The data type of the input (fasta, peptide list)
  -l, --length {8,9,10,11,12,13,14,15,16,17}
                        The length of peptides
  -a ALLELES, --alleles ALLELES
                        Path to the allele file (one per line in new
                        nomenclature)
  -op OPTIONS, --options OPTIONS
                        Additional options that get directly past to the tool
  -o OUTPUT, --output OUTPUT
                        Path to the output file

"""
import sys
import pandas
import argparse

from Fred2.Core import Protein, Peptide, Allele
from Fred2.IO import read_fasta
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.Core import generate_peptides_from_proteins

def read_lines(file, in_type=Peptide):
    peptides = []

    with open(file, "r") as f:
        for l in f:
            if not l.startswith("#") and l.strip() != "" and not l.startswith("Epitope") and not l.startswith("Sequence"):
                pep = l.split()[0].strip()
                peptides.append(in_type(pep))
    return peptides

def main():
    #Specify CTD interface
    # Every CTD Model has to have at least a name and a version, plus any of the optional attributes below them.
    model = argparse.ArgumentParser(description='Process some integers.')

    model.add_argument('-m',
        '--method',
        type=str,
        choices=EpitopePredictorFactory.available_methods().keys(),
        default="bimas",
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
        choices=["fasta","peptide"],
        type=str,
        default="fasta",
        help='The data type of the input (fasta, peptide list)'
        )

    model.add_argument('-l',
        '--length',
        choices=range(8, 18),
        type=int,
        default=9,
        help='The length of peptides'
        )

    model.add_argument('-a',
        '--alleles',
        type=str,
        required=True,
        help='Path to the allele file (one per line in new nomenclature)'
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
        peptides = generate_peptides_from_proteins(proteins, args.length)
    elif args.type == "peptide":
        peptides = read_lines(args.input, in_type=Peptide)
    else:
        sys.stderr.write('Input type not known\n')
        return -1

    #read in alleles
    alleles = read_lines(args.alleles, in_type=Allele)
    if args.version == "":
        result = EpitopePredictorFactory(args.method).predict(peptides, alleles, options=args.options)
    else:
        result = EpitopePredictorFactory(args.method, version=args.version).predict(peptides, alleles,
                                                                 options=args.options)

    #write to TSV columns sequence method allele-scores...,protein-id/transcript-id
    with open(args.output, "w") as f:
        proteins = "\tAntigen ID" if args.type == "fasta" else ""
        alleles = result.columns
        f.write("Sequence\tMethod\t"+"\t".join(a.name for a in alleles)+proteins+"\n")
        for index, row in result.iterrows():
            p = index[0]
            method = index[1]
            proteins =  "\t"+",".join( prot.transcript_id for prot in p.get_all_proteins()) if args.type == "fasta" else ""
            f.write(str(p)+"\t"+method+"\t"+"\t".join("%.3f"%row[a] for a in alleles)+proteins+"\n")

    return 0

if __name__ == "__main__":
    sys.exit(main())
