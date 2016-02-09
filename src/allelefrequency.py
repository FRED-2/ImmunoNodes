#!/usr/bin/env python
"""
Command line tool to assign allele frequencies to alleles lists

usage: allelefrequency.py [-h]
                          [-g {Oceania,South-East_Asia,South-West_Asia,Europe,North_Africa,Australia,South_America,
                               Sub-Saharan_Africa,Other,North-East_Asia,North_America}]
                          [-p {Moroccan_99,Moroccan_98,Muong,Singapore_Chinese),Finn_90,North_America_(Eu),Tamil,
                               Israeli_Jews,Filipino,Omani,Turk,Cuban_(Eu),Hakka,Zambian,Kenyan_142,North_America_(As),
                               Yuendumu,Algerian_99,Slovenian,Saisiat,North_America_(Hi),Rwandan,Okinawan,Mandenka,Zulu,
                               Georgian,Ivatan,Ugandan,Thai,Han-Chinese_572,Thao,Sioux,Bulgarian,Kinh,South_Indian,
                               Atayal,Kenyan_Lowlander,Metalsa,Chaouya,PNG_Highlander,Chinese,Puyuma_49,Tsou,Bari,
                               Buriat,Minnan,Cuban_(Af_Eu),Brazilian,Tuva,Canoncito,Kenyan_Highlander,East_Timorese,
                               Malay,Zuni,PNG_Lowlander_48,North_America_(Af),Ticuna,New_Delhi,Bunun,Brazilian_(Af_Eu,
                               Pima_99,Ami_97,Korean_200,Irish,Czech,Guarani-Kaiowa,Toroko,Ryukuan,Pima_17,
                               Arab_Druze,Moluccan,Shona,Yupik,Javanese_Indonesian,Yami,Cape_York,Maya,Paiwan_51,
                               Guarani-Nandewa,Rukai,Kimberley,American_Samoa,Kurdish,Siraya,Groote_Eylandt,
                               Han-Chinese_149,Central_America,Pazeh,Seri,Mexican,Lacandon,PNG_Lowlander_95,
                               Amerindian,Croatian,Doggon}]
                          [-a ALLELES] [-t THRESHOLD] -o OUTPUT

Commandline tool for allele frequency assignment

optional arguments:
  -h, --help            show this help message and exit
  -g  GEOGRAPHIC --geographic {Oceania,South-East_Asia,South-West_Asia,Europe,North_Africa,Australia,South_America,
                               Sub-Saharan_Africa,Other,North-East_Asia,North_America}
                        Specifies the geographic region (e.g. Europe)
  -p POPULATION, --population {Moroccan_99,Moroccan_98,Muong,Singapore_(Chinese),Finn_90,North_America_(Eu),Tamil,
                               Israeli_Jews,Filipino,Omani,Turk,Cuban_(Eu),Hakka,Zambian,Kenyan_142,North_America_(As),
                               Yuendumu,Algerian_99,Slovenian,Saisiat,North_America_(Hi),Rwandan,Okinawan,Mandenka,Zulu,
                               Georgian,Ivatan,Ugandan,Thai,Han-Chinese_572,Thao,Sioux,Bulgarian,Kinh,South_Indian,
                               Atayal,Kenyan_Lowlander,Metalsa,Chaouya,PNG_Highlander,Chinese,Puyuma_49,Tsou,Bari,
                               Buriat,Minnan,Cuban_(Af_Eu),Brazilian,Tuva,Canoncito,Kenyan_Highlander,East_Timorese,
                               Malay,Zuni,PNG_Lowlander_48,North_America_(Af),Ticuna,New_Delhi,Bunun,Brazilian_(Af_Eu),
                               Pima_99,Ami_97,Korean_200,Irish,Czech,Guarani-Kaiowa,Toroko,Ryukuan,Pima_17,Arab_Druze,
                               Moluccan,Shona,Yupik,Javanese_Indonesian,Yami,Cape_York,Maya,Paiwan_51,Guarani-Nandewa,
                               Rukai,Kimberley,American_Samoa,Kurdish,Siraya,Groote_Eylandt,Han-Chinese_149,
                               Central_America,Pazeh,Seri,Mexican,Lacandon,PNG_Lowlander_95,Amerindian,Croatian,Doggon}
                        Specifies the population (e.g. Irish)
  -a ALLELES, --alleles ALLELES
                        Path to the allele file (one per line in new
                        nomenclature)
  -t THRESHOLD, --threshold THRESHOLD
                        Frequency threshold to include alleles
  -o OUTPUT, --output OUTPUT
                        Path to the output file


"""
from Fred2.Core import Allele
from Fred2.IO import read_lines

import sys
from data.geo import geo
from data.pop import pop
import argparse


def main():
    model = argparse.ArgumentParser(
        description='Commandline tool for allele frequency assignment',
        )

    model.add_argument(
        '-g','--geographic',
        choices=geo.keys(),
        type=str,
        default="",
        help='Specifies the geographic region (e.g. Europe)',
        )
    model.add_argument(
        '-p','--population',
        choices=pop.keys(),
        type=str,
        default="",
        help='Specifies the population (e.g. Irish)',
        )
    model.add_argument(
        '-a','--alleles',
        type=str,
        default="",
        help='Path to the allele file (one per line in new nomenclature)',
        )

    model.add_argument(
        '-t','--threshold',
        type=float,
        default=0.001,
        help="Frequency threshold to include alleles",
    )

    model.add_argument(
        '-o','--output',
        type=str,
        required=True,
        help='Path to the output file',
        )

    args = model.parse_args()

    if (args.population != "" and args.geographic != "") or \
            (args.population == "" and args.geographic == ""):
        sys.stderr.write("Please select either a geographic region or a population.\n")
        return -1

    thr = float(args.threshold)
    with open(args.output, "w") as f:
        if args.geographic in geo:
            freqs = geo[args.geographic]
        elif args.population in pop:
            freqs = pop[args.population]
        else:
            sys.stderr.write("{pop} could not be found\n".format(pop=args.population))
            return -1

        if args.alleles != "":
            alleles = read_lines(args.alleles, in_type=Allele)
            for a in alleles:
                fr = freqs.get(a.name, 0.0)
                if fr >= thr:
                    f.write(a.name+"\t"+"%.3f"%fr+"\n")
        else:
            for a, fr in freqs.iteritems():
                fr = fr
                if fr >= thr:
                    f.write(a+"\t"+"%.3f"%fr+"\n")
    return 0

if __name__ == "__main__":
    sys.exit(main())