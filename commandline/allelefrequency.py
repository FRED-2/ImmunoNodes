"""
Command line tool to assign allele frequencies to alleles lists

"""
from CTDopts.CTDopts import CTDModel
from Fred2.Core import Allele
from Fred2.IO import read_lines

import sys
from data.geo import geo
from data.pop import pop


def main():
    #Specify CTD interface
    # Every CTD Model has to have at least a name and a version, plus any of the optional attributes below them.
    model = CTDModel(
        name='AlleleFrequency',  # required
        version='1.0',  # required
        description='Commandline tool for allele frequency assignment',
        manual='manual string',
        executableName='allelefrequency',
        )

    model.add(
        'population',
        choices=geo.keys()+pop.keys(),
        type=str,
        required=True,
        description='Specifies the population (e.g. Europe)',
        short_name="p"
        )

    model.add(
        'alleles',
        type="input-file",
        default=None,
        description='Path to the allele file (one per line in new nomenclature)',
        short_name="a"
        )

    model.add(
        'threshold',
        type=float,
        num_range=(0.0,1.0),
        default=0.001,
        description="Frequency threshold to include alleles",
        short_name="t"
    )

    model.add(
        'output',
        type="output-file",
        required=True,
        description='Path to the output file',
        short_name="o"
        )

    model.add(
        'ctdout',
        default=None,
        type="output-file",
        description='Output path to for cds',
        short_name="cds"
        )

    args_str = sys.argv[1:] if sys.argv[1:] else ["--help"]
    args = model.parse_cl_args(cl_args=args_str)
    print args
    if args["ctdout"] is not None:
        model.write_ctd(args[args["ctdout"]])
        return 0

    thr = float(args["threshold"])
    with open(args["output"], "w") as f:
        if args["population"] in geo:
            freqs = geo[args["population"]]
        elif args["population"] in pop:
            freqs = geo[args["population"]]
        else:
            sys.stderr("{pop} could not be found".format(pop=args["population"]))
            return -1

        if args["alleles"] is not None:
            alleles = read_lines(args["alleles"], type=Allele)
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