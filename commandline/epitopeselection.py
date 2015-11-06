#!/usr/bin/env python
"""
Command line tool for epitope selection
"""
import sys
import pandas

from CTDopts.CTDopts import CTDModel

from Fred2.EpitopeSelection.OptiTope import OptiTope
from Fred2.Core import Allele, Peptide, Protein


def generate_epitope_result(input, alleles):
    """
    generates EpitopePredictionResult from output of epitopeprediction and neoepitopeprediction
    """
    #first generate alleles in allele file
    alleles = {}
    with open(alleles, "r") as af:
        for l in af:
            allele, freq = l.split("\t")
            alleles[allele] = Allele(allele, prob=float(freq))

    r_raw = pandas.read_csv(input, sep="\t")
    print r_raw
    sys.exit()


def to_csv(out_file, result, instance, pred_method):
    """
        Writes model to CSV
    """
    with open(out_file, "w") as f:
        f.write("Prediction method: " + pred_method + "\n\n")

        cons = ["Maximum number of epitopes to select = " + str(int(instance.k)) + "\n"]
        if float(instance.t_c) > 0:
            cons.append("Epitope conservation =>" + str(float(instance.t_c) * 100) + "%\n")

        if float(instance.t_allele) > 0:
            cons.append("Covered alleles >= " + str(int(instance.t_allele)) + "\n")

        if float(instance.t_var) > 0:
            cons.append("Covered antigens >= " + str(int(instance.t_var)) + "\n")
        f.write("CONSTRAINTS\n" + "".join(cons) + "\n")

        res = ["Selected epitopes," + str(len(result)) + ""]

        if int(instance.t_var) > 0:
            cov_anti = []
            for an in instance.Q:
                for e in result:
                    if e in instance.E_var[an]:
                        cov_anti.append(an)
            cov_anti = set(cov_anti)
            res.append("Covered antigens," + str(len(cov_anti)) + " of " + str(len(instance.Q)) + "")
        cov_als = []
        res_set = set(result)
        locus = {}
        for a in instance.A:
            eps_of_all_i = list(instance.A_I[a])
            if res_set.intersection(set(eps_of_all_i)):
                cov_als.append(a)
                locus.setdefault(str(a).split("*")[0], set()).add(a)
        cov_als = set(cov_als)
        res.append("Covered alleles," + str(len(cov_als)) + " of " + str(len(instance.A)) + "")
        res.append("Locus coverage:\n")

        pop_cov = 1
        for k, g in locus.iteritems():
            locus = list(g)
            pop_cov *= (1.0 - sum(float(instance.p[a]) for a in locus)) ** 2
            covered = len(locus) / sum(1 for a in instance.A if a.split("*")[0] == k)
            res.append(",%s,%.2f" % (k, covered * 100))
        res.append("Population converage:,,%.2f" % ((1.0 - pop_cov) * 100))
        f.write("RESULTS\n" + "\n".join(res) + "\n\n")

        is_antigen_cons = int(instance.t_var) > 0
        header = "Epitope,Conservation,Fraction of overall immunogenicity,Covered alleles%s\n" % (
            ",Covered antigens" if is_antigen_cons else "" )

        rows = []
        overall_imm = sum(float(instance.i[e, a]) * float(instance.p[a]) for e in result for a in instance.A)
        for e in result:
            row = str(e) + ","
            if float(instance.t_c) > 0:
                row += str(float(instance.c[e]) * 100) + ","
            else:
                row += "100%,"
            row += "%0.2f," % (sum(float(instance.i[e, a]) * float(instance.p[a]) for a in instance.A) / overall_imm)
            row += "%s" % " ".join(str(a) for a in instance.A if e in instance.A_I[a])
            # row += "%s"%"\n".join( str(float(instance.i[e,a])*float(instance.p[a]))  for a in instance.A if e in instance.A_I[a] )
            if is_antigen_cons:
                row += ",%s" % " ".join(str(q) for q in instance.Q if e in instance.E_var[q])
            # row+="\n"
            rows.append(row)
        f.write(header + "\n".join(rows) + "\n\n")
        f.write("TARGET POPULATION/INDIVIDUAL\n")
        f.write(",Allele,Probability\n")
        for a in instance.A:
            f.write(",%s,%.5f\n" % (str(a), float(instance.p[a])))


def main():
    '''
        some input stuff
    '''
    parser = CTDModel(name="EpitopeSelection",
                      version='1.0',  # required
                      description="Epitope Selection for vaccine design.",
                      manual='manual string',
                      executableName='epitopeselection'
    )
    parser.add("input",
               required=True,
               type="input-file",
               help="Peptide with immunogenicity file (from epitopeprediction)",
               short_name="i"
    )
    parser.add("alleles",
               required=True,
               type="input-file",
               help="Allele file with frequencies (one allele and frequency per line)",
               short_name="a"
    )

    parser.add("k",
               required=False,
               type=int,
               default=10,
               help="Specifies the number of epitopes to select",
               short_name="k")
    parser.add("threshold",
               required=True,
               type=float,
               default=0.,
               help="Specifies the binding threshold for all alleles",
               short_name="thr")
    parser.add("output",
               required=True,
               type="output-file",
               help="Specifies the output path. Results will be written to CSV",
               short_name="o"
    )

    parser.add("cons-allele",
               required=False,
               type=float,
               default=None,
               help="Activates allele coverage constraint with specified threshold",
               short_name="cal")
    parser.add("cons-antigen",
               required=False,
               type=float,
               default=None,
               help="Activates antigen coverage constraint with specified threshold",
               short_name="can")
    parser.add("cons-conservation",
               required=False,
               type=float,
               help="Activates conservation constraint with specified threshold",
               default=None,
               short_name="ccon")
    parser.add("conservation",
               default=None,
               required=False,
               type="input-file",
               help="Specifies a Conservation file. First column is the peptide seq second column the conservation.",
               short_name="cons")
    parser.add(
        'ctdout',
        default=None,
        type="output-file",
        description='Output path to for cds',
        short_name="ctd"
    )

    args_str = sys.argv[1:] if sys.argv[1:] else ["--help"]
    args = parser.parse_cl_args(cl_args=args_str)
    print args
    if args["ctdout"] is not None:
        parser.write_ctd(args[args["ctdout"]])
        return 0

    epitopePrediciton, method = generate_epitope_result(args["input"], args["alleles"])
    thresh = {a.name: args["threshold"] for a in epitopePrediciton.columns}
    opti = OptiTope(epitopePrediciton, threshold=thresh, k=args.k, solver="cbc", threads=1, verbosity=0)

    # set constraints
    if args["cons-allele"] > 0.0:
        #print "allele constraiont enforced"
        opti.activate_allele_coverage_const(args.cons_allele / 100)

    if args["cons-antigen"] > 0.0:
        opti.activate_antigen_coverage_const(args.cons_antigen / 100)

    if args["cons-conservation"] > 0.0:

        if args["conservation"] is not None:
            conservation = {}
            with open(args["conservation"], "r") as f:
                for l in f:
                    seq, cons = l.replace(",", " ").replace(";", " ").split()
                    conservation[seq.strip().upper()] = float(cons.strip())
            opti.activate_epitope_conservation_const(args["cons-conservation"] / 100.0, conservation=conservation)
        else:
            opti.activate_epitope_conservation_const(args["cons-conservation"] / 100.0)
    try:
        result = opti.solve()
        #print result
        to_csv(args["output"], result, opti.instance, method)
        return 0
    except ValueError as e:
        sys.stderr.write("Could not optimally solve the problem. Please modify your constraints.\n"+str(e))
        return -1
    except Exception as e:
        sys.stderr.write(str(e))
        return -1


if __name__ == "__main__":
    sys.exit(main())