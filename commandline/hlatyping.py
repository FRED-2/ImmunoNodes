"""
Commandline tool for hla typing prediction


"""
import sys
import tempfile

from CTDopts import CTDModel

from Fred2.HLAtyping import HLATypingFactory


def main():
#Specify CTD interface
    # Every CTD Model has to have at least a name and a version, plus any of the optional attributes below them.
    model = CTDModel(
        name='EpitopePredicton',  # required
        version='1.0',  # required
        description='Commandline tool for epitope prediction',
        manual='manual string',
        executableName='epitopeprediction',
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
        type="input_file",
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
        type="output_file",
        description='Path to the output file'
        )

    args = model.parse_cl_args(cl_args=sys.argv[1:])

    #fasta protein
    genotype = HLATypingFactory(args["name"], version=args["version"]).predict(args["input"], args["tmp_output"],
                                                                               options=args["options"])
    with open(args["output"], "w") as f:
        f.write("\n".join(a.name for a in genotype))

    return 0


if __name__ == "__main__":
    sys.exit(main())