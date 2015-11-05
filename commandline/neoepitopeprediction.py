#!/usr/bin/env python
"""
Generator for polymorphic epitopes form annovar exonic files command line tool

"""
import sys
import tempfile
import subprocess
import os

from CTDopts.CTDopts import CTDModel

from Fred2.Core import Allele, generate_peptides_from_variants
from Fred2.IO import read_lines, read_annovar_exonic
from Fred2.IO.MartsAdapter import MartsAdapter
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.Core.Variant import VariationType

ANOVA_PREP = "convert2annovar.pl -format vcf4 -withzyg %s > %s "
ANOVA = "annotate_variation.pl -out %s -build %s %s /abi-projects/etk/software/annovar/humandb/ -dbtype ensGene"

#ensGene, refGene
#TODO: ANNOVAR PATHS and additional files have to be generic
#TODO: MartsAadapter has to changed underlying REST server depending on reference used
#TODO: how should the output format look like?

def main():
    #Specify CTD interface
    # Every CTD Model has to have at least a name and a version, plus any of the optional attributes below them.
    model = CTDModel(
        name='Polymorphic Epitope Prediction',  # required
        version='1.0',  # required
        description='Commandline tool for Neo-epitope prediction',
        manual='manual string',
        executableName='neoepitopeprediction',
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
        description='Path to the input file (vcf format)'
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

    model.add("snps",
              required=False,
              type=bool,
              default=True,
              help="Filter for variations (excluding snps)"
    )

    model.add("indels",
              required=False,
              type=bool,
              default=False,
              help="Filter for variations (excluding indels)"
    )

    model.add("frame-shift",
              required=False,
              type=bool,
              default=False,
              help="Filter for variations (excluding fram-shifts)"
    )


    model.add("reference",
              required=False,
              type=str,
              choices=["hg19"],
              default="hg19",
              help="Human reference for ANOVAR annotation (default: hg19)"
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
    print args
    if args["ctdout"] is not None:
        model.write_ctd(args[args["ctdout"]])
        return 0

    #generate polymorphic peptides with annovar annotation
    tmp_annovar_input = tempfile.NamedTemporaryFile(delete=False)
    tmp_annovar_input.close()
    subprocess.call(ANOVA_PREP%(args["input"], tmp_annovar_input.name), shell=True)

    #call anova main annotation function with prepared input file
    #print ANOVA%(tmp_annovar_input.name,tmp_annovar_input.name)
    subprocess.call(ANOVA%(args["reference"], tmp_annovar_input.name, tmp_annovar_input.name), shell=True)
    anovar_exonic = tmp_annovar_input.name+".exonic_variant_function"
    os.remove(tmp_annovar_input.name+".variant_function")
    os.remove(tmp_annovar_input.name)


    vars = filter(lambda x: x.type != VariationType.UNKNOWN, read_annovar_exonic(anovar_exonic))
    os.remove(anovar_exonic)
    #here we could filter the vars
    #print "Variants", vars

    if not args["snps"]:
        vars = filter(lambda x: x.type != VariationType.SNP, vars)

    if not args["indels"]:
        vars = filter(lambda x: x.type not in [VariationType.INS, VariationType.DEL,
                                               VariationType.FSINS, VariationType.FSDEL], vars)

    if not args["frame-shift"]:
        vars = filter(lambda x: x.type not in [VariationType.FSINS, VariationType.FSDEL], vars)

    if not vars:
        #to_html_error(args.output, "No variants could be found")
        sys.stderr.write('No variants could be found. Please check whether your input is of human origin.\n')
        return -1

    #TODO: filter for specific genes of interrest -> should be done in generate_peptides_from_vars
    #->less work!
    peptides = generate_peptides_from_variants(vars, args.length, MartsAdapter())

    #read in alleles
    alleles = read_lines(args["alleles"], type=Allele)
    if args["version"] is None:
        result = EpitopePredictorFactory(args["method"]).predict(peptides, alleles, options=args["options"])
    else:
        result = EpitopePredictorFactory(args["method"], version=args["version"]).predict(peptides, alleles,
                                                                 options=args["options"])

    #write to TSV columns sequence method allele-scores...,protein-id/transcript-id
    with open(args["output"], "w") as f:
        proteins = "\tProtein ID" if args["type"] == "fasta" else ""
        alleles = result.columns
        f.write("Sequence\tMethod\t"+"\t".join(a.name for a in alleles)+proteins+"\n")
        for index, row in result.iterrows():
            p = index[0]
            method = index[1]
            proteins = ",".join( prot.transcript_id for prot in p.get_all_proteins()) if args["type"] == "fasta" else ""
            f.write(str(p)+"\t"+method+"\t"+"\t".join("%.3f"%row[a] for a in alleles)+"\t"+proteins+"\n")
    return 0

if __name__ == "__main__":
    sys.exit(main())