#!/usr/bin/env python
"""
Commandline tool for hla typing prediction


"""
import sys
import tempfile

from CTDopts.CTDopts import CTDModel

from Fred2.HLAtyping import HLATypingFactory


def main():
#Specify CTD interface
    # Every CTD Model has to have at least a name and a version, plus any of the optional attributes below them.
    model = CTDModel(
        name='HLATyping',  # required
        version='1.0',  # required
        description='Commandline tool for HLA typing',
        manual='manual string',
        executableName='hlatyping',
        )

    model.add(
        'method',
        type=str,
        choices=HLATypingFactory.available_methods().keys(),
        default="optitype",
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
        'options',
        type=str,
        default=None,
        description="Additional options that get directly past to the tool"
    )

    model.add(
        'tmp_output',
        type="str",
        default=tempfile.tempdir,
        description='Path to a temporary output file or directory used by the HLAtyping tool'
        )

    model.add(
        'output',
        type="output-file",
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
    genotype = HLATypingFactory(args["name"], version=args["version"]).predict(args["input"], args["tmp_output"],
                                                                               options=args["options"])
    with open(args["output"], "w") as f:
        f.write("\n".join("HLA-"+a.name for a in genotype))

    return 0


if __name__ == "__main__":
    sys.exit(main())