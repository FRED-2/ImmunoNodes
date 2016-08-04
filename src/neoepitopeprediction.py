#!/usr/bin/env python
"""
Commandline tool for (neo)epitope prediction

usage: neoepitopeprediction.py [-h]
                               [-m {netmhc,smmpmbec,syfpeithi,netmhcpan,netctlpan,smm,tepitopepan,arb,pickpocket,epidemix,netmhcii,netmhciipan,comblibsidney,unitope,hammer,svmhc,bimas}]
                               [-v VCF] [-t {VEP,ANNOVAR}] [-p PROTEINS]
                               [-l {8,9,10,11,12,13,14,15,16,17}] -a ALLELES
                               [-r REFERENCE] [-fINDEL] [-fFS] [-fSNP] -o
                               OUTPUT

Neoepitope prediction for TargetInsepctor.

optional arguments:
  -h, --help            show this help message and exit
  -m {netmhc,smmpmbec,syfpeithi,netmhcpan,netctlpan,smm,tepitopepan,arb,pickpocket,epidemix,netmhcii,netmhciipan,comblibsidney,unitope,hammer,svmhc,bimas}, --method {netmhc,smmpmbec,syfpeithi,netmhcpan,netctlpan,smm,tepitopepan,arb,pickpocket,epidemix,netmhcii,netmhciipan,comblibsidney,unitope,hammer,svmhc,bimas}
                        The name of the prediction method
  -v VCF, --vcf VCF     Path to the vcf input file
  -t {VEP,ANNOVAR}, --type {VEP,ANNOVAR}
                        Type of annotation tool used (Variant Effect
                        Predictor, ANNOVAR exonic gene annotation)
  -p PROTEINS, --proteins PROTEINS
                        Path to the protein ID input file (in HGNC-ID)
  -l {8,9,10,11,12,13,14,15,16,17}, --length {8,9,10,11,12,13,14,15,16,17}
                        The length of peptides
  -a ALLELES, --alleles ALLELES
                        Path to the allele file (one per line in new
                        nomenclature)
  -r REFERENCE, --reference REFERENCE
                        The reference genome used for varinat annotation and
                        calling.
  -fINDEL, --filterINDEL
                        Filter insertions and deletions (including
                        frameshifts)
  -fFS, --filterFSINDEL
                        Filter frameshift INDELs
  -fSNP, --filterSNP    Filter SNPs
  -o OUTPUT, --output OUTPUT
                        Path to the output file

"""
import time
import sys
import argparse

from Fred2.Core import Protein, Peptide, Allele, MutationSyntax, Variant
from Fred2.Core.Variant import VariationType
from Fred2.IO import read_lines, MartsAdapter, read_annovar_exonic
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.Core import generate_peptides_from_proteins, generate_peptides_from_variants
from Fred2.IO.ADBAdapter import EIdentifierTypes, EAdapterFields


MARTDBURL = {"GRCH37": "http://grch37.ensembl.org/biomart/martservice?query=",
             "GRCH38": "http://www.ensembl.org/biomart/martservice?query="} # is corrently set to GRCh38

def read_variant_effect_predictor(file, gene_filter=None):
    """
    Reads a VCF (v4.1) file generatede by variant effect predictor and generates variant objects

    :param str file: Path to vcf file
    :param list gene_filter: List of proteins (in HGNC) of inerrest. Variants are filter according to this list
    :return: list(Variant) - a list of Fred2.Core.Variant objects
    """
    vars = []
    def get_type(ref, alt):
        """
            returns the variant type
        """
        if len(ref)==1 and len(alt)==1:
            return VariationType.SNP
        if len(ref)>0 and len(alt)==0:
            if len(ref)%3 == 0:
                return VariationType.DEL
            else:
                return VariationType.FSDEL
        if len(ref) == 0 and len(alt)>0:
            if len(alt)% 3 == 0:
                return VariationType.INS
            else:
                return VariationType.FSINS
        return VariationType.UNKNOWN

    coding_types = set(["3_prime_UTR_variant", "5_prime_UTR_variant", "start_lost", "stop_gained",
        "frameshift_variant", "start_lost", "inframe_insertion", "inframe_deletion", "missense_variant",
        "protein_altering_variant", "splice_region_variant", "incomplete_terminal_codon_variant", "stop_retained_variant",
        "synonymous_variant", "coding_sequence_variant"])

    with open(file, "r") as f:
        for i,l in enumerate(f):

            #skip comments
            if l.startswith("#") or l.strip() == "":
                continue

            chrom, gene_pos,var_id,ref,alt,_,filter_flag,info= l.strip().split("\t")[:8]
            coding = {}
            isSynonymous = False

            for co in info.split(","):
                #Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|SYMBOL|SYMBOL_SOURCE|HGNC_ID|CCDS">
                _,gene,transcript_id,transcript_type,var_type,_,transcript_pos,prot_pos,_,_,_,distance,strand,HGNC_ID = co.strip().split("|")[:14]

                #pass every other feature type except Transcript (RegulatoryFeature, MotifFeature.)
                #pass genes that are uninterresting for us
                if transcript_type != "Transcript" or (HGNC_ID not in gene_filter and gene_filter):
                    continue

                #pass all intronic and other mutations that do not directly influence the protein sequence
                if any(t in coding_types for t in var_type.split("&")):
                    #generate mutation syntax

                    #positioning in Fred2 is 0-based!!!
                    if transcript_pos != "":
                        coding[transcript_id] = MutationSyntax(transcript_id, int(transcript_pos.split("-")[0])-1, 
                            -1 if prot_pos  == "" else int(prot_pos.split("-")[0])-1, co, "", geneID=HGNC_ID)

                #is variant synonymous?
                isSynonymous = any(t == "synonymous_variant" for t in var_type.split("&"))
            if coding:
                vars.append(Variant(var_id, get_type(ref, alt), chrom, int(gene_pos), ref.upper(), alt.upper(), coding, False, isSynonymous))
    return vars


def main():

    model = argparse.ArgumentParser(description='Neoepitope prediction for TargetInsepctor.')

    model.add_argument(
        '-m','--method',
        type=str,
        choices=EpitopePredictorFactory.available_methods().keys(),
        default="bimas",
        help='The name of the prediction method'
        )


    model.add_argument(
        '-v', '--vcf',
        type=str,
        default=None,
        help='Path to the vcf input file'
        )

    model.add_argument(
        '-t', '--type',
        type=str,
        choices=["VEP", "ANNOVAR"],
        default="VEP",
        help='Type of annotation tool used (Variant Effect Predictor, ANNOVAR exonic gene annotation)'
        )

    model.add_argument(
        '-p','--proteins',
        type=str,
        default=None,
        help='Path to the protein ID input file (in HGNC-ID)'
        )

    model.add_argument(
        '-l','--length',
        choices=range(8, 18),
        type=int,
        default=9,
        help='The length of peptides'
        )

    model.add_argument(
        '-a','--alleles',
        type=str,
        required=True,
        help='Path to the allele file (one per line in new nomenclature)'
        )

    model.add_argument(
        '-r' ,'--reference',
        type=str,
        default='GRCh38',
        help='The reference genome used for varinat annotation and calling.'
        )

    model.add_argument(
        '-fINDEL' ,'--filterINDEL',
        action="store_true",
        help='Filter insertions and deletions (including frameshifts)'
        )

    model.add_argument(
        '-fFS' ,'--filterFSINDEL',
        action="store_true",
        help='Filter frameshift INDELs'
        )

    model.add_argument(
        '-fSNP' ,'--filterSNP',
        action="store_true",
        help='Filter SNPs'
        )

    model.add_argument(
        '-o','--output',
        type=str,
        required=True,
        help='Path to the output file'
        )
    model.add_argument(
        '-etk','--etk',
        action="store_true",
        help=argparse.SUPPRESS
        )

    args = model.parse_args()

    martDB = MartsAdapter(biomart=MARTDBURL[args.reference.upper()])
    transcript_to_genes = {}

    if args.vcf is None and args.proteins is None:
        sys.stderr.write("At least a vcf file or a protein id file has to be provided.\n")
        return -1

    # if vcf file is given: generate variants and filter them if HGNC IDs ar given
    if args.vcf is not None:
        protein_ids = []
        if args.proteins is not None:
            with open(args.proteins, "r") as f:
                for l in f:
                    l = l.strip()
                    if l != "":
                        protein_ids.append(l)
        if args.type == "VEP":
            print "in here"
            variants = read_variant_effect_predictor(args.vcf, gene_filter=protein_ids)
            print variants
        else:
            variants = read_annovar_exonic(args.vcf, gene_filter=protein_ids)

        variants = filter(lambda x: x.type != VariationType.UNKNOWN, variants)

        if args.filterSNP:
            variants = filter(lambda x: x.type != VariationType.SNP, variants)

        if args.filterINDEL:
            variants = filter(lambda x: x.type not in [VariationType.INS,
                                                       VariationType.DEL,
                                                       VariationType.FSDEL,
                                                       VariationType.FSINS], variants)

        if args.filterFSINDEL:
            variants = filter(lambda x: x.type not in [VariationType.FSDEL, VariationType.FSINS], variants)

        if not variants:
            sys.stderr.write("No variants left after filtering. Please refine your filtering criteria.\n")
            return -1

        epitopes = filter(lambda x:any(x.get_variants_by_protein(tid) for tid in x.proteins.iterkeys()),
                        generate_peptides_from_variants(variants,
                                                int(args.length), martDB, EIdentifierTypes.ENSEMBL))

        for v in variants:
            for trans_id,coding in v.coding.iteritems():
                transcript_to_genes[trans_id] = coding.geneID

    #else: generate protein sequences from given HGNC IDs and than epitopes
    else:
        proteins = []
        with open(args.proteins, "r") as f:
            for l in f:
                ensembl_ids = martDB.get_ensembl_ids_from_id(l.strip(), type=EIdentifierTypes.HGNC)[0]
                protein_seq = martDB.get_product_sequence(ensembl_ids[EAdapterFields.PROTID])
                if protein_seq is not None:
                    transcript_to_genes[ensembl_ids[EAdapterFields.TRANSID]] = l.strip()
                    proteins.append(Protein(protein_seq, gene_id=l.strip(), transcript_id=ensembl_ids[EAdapterFields.TRANSID]))
        epitopes = generate_peptides_from_proteins(proteins, int(args.length))


    #read in allele list
    alleles = read_lines(args.alleles, in_type=Allele)

    result = EpitopePredictorFactory(args.method).predict(epitopes, alleles=alleles)

    with open(args.output, "w") as f:
        alleles = result.columns
        var_column = " Variants" if args.vcf is not None else ""
        f.write("Sequence\tMethod\t"+"\t".join(a.name for a in alleles)+"\tAntigen ID\t"+var_column+"\n")
        for index, row in result.iterrows():
            p = index[0]
            method = index[1]
            proteins = ",".join(set([transcript_to_genes[prot.transcript_id.split(":FRED2")[0]] for prot in p.get_all_proteins()]))
            vars_str = ""

            if args.vcf is not None:
                vars_str = "\t"+"|".join(set(prot_id.split(":FRED2")[0]+":"+",".join(repr(v) for v in set(p.get_variants_by_protein(prot_id)))
                                                                            for prot_id in p.proteins.iterkeys()
                                          if p.get_variants_by_protein(prot_id)))
            
            f.write(str(p)+"\t"+method+"\t"+"\t".join("%.3f"%row[a] for a in alleles)+"\t"+proteins+vars_str+"\n")

    if args.etk:
        with open(args.output.rsplit(".",1)[0]+"_etk.tsv", "w") as g:
            alleles = result.columns
            g.write("Alleles:\t"+"\t".join(a.name for a in alleles)+"\n")
            for index, row in result.iterrows():
                p = index[0]
                proteins = " ".join(set([transcript_to_genes[prot.transcript_id.split(":FRED2")[0]] for prot in p.get_all_proteins()]))
                g.write(str(p)+"\t"+"\t".join("%.3f"%row[a] for a in alleles)+"\t"+proteins+"\n")
    return 0



if __name__ == "__main__":
    sys.exit(main())

