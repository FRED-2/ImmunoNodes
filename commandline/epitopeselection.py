#!/usr/bin/env python
"""
Command line tool for epitope selection
"""
import sys
import pandas
import collections
from CTDopts.CTDopts import CTDModel

from Fred2.EpitopeSelection.OptiTope import OptiTope
from Fred2.Core import Allele, Peptide, Protein, EpitopePredictionResult


def generate_epitope_result(input, allele_file):
    """
    generates EpitopePredictionResult from output of epitopeprediction and neoepitopeprediction
    """
    #first generate alleles in allele file
    alleles = {}
    with open(allele_file, "r") as af:
        for l in af:
            allele, freq = l.split("\t")
            alleles[allele] = Allele(allele, prob=float(freq))

    r_raw = pandas.read_csv(input, sep="\t")
    res_dic = {}
    method = r_raw.loc[0, "Method"]
    columns = set(["Sequence", "Method", "Protein ID", "Variant"])
    alleles_raw = [c for c in r_raw.columns if c not in columns]
    for k, row in r_raw.iterrows():
        seq = row["Sequence"]
        protPos = collections.defaultdict(list)
        try:
            protPos = {Protein("", _gene_id=p, _transcript_id=p): [0] for p in str(row["Protein ID"]).split(",")}
        except KeyError:
            pass
        for a in alleles_raw:
            if a in alleles:
                if alleles[a] not in res_dic:
                    res_dic[alleles[a]] = {}
                res_dic[alleles[a]][Peptide(seq, protein_pos=protPos)] = float(row[a])

    if not res_dic:
        sys.stderr.write("HLA alleles of population and HLA used for prediction did not overlap.")
        sys.exit(-1)

    df_result = EpitopePredictionResult.from_dict(res_dic)
    df_result.index = pandas.MultiIndex.from_tuples([tuple((i, method)) for i in df_result.index],
                                                        names=['Seq', 'Method'])
    return df_result, method


def to_csv(out_file, result, instance, pred_method):
    """
        Writes model to CSV
    """
    with open(out_file, "w") as f:
        f.write("Prediction method: " + pred_method + "\n\n")

        cons = ["Maximum number of epitopes to select = " + str(int(instance.k.value)) + "\n"]
        if float(instance.t_c.value) > 0:
            cons.append("Epitope conservation =>" + str(float(instance.t_c.value) * 100) + "%\n")

        if float(instance.t_allele.value) > 0:
            cons.append("Covered alleles >= " + str(int(instance.t_allele.value)) + "\n")

        if float(instance.t_var.value) > 0:
            cons.append("Covered antigens >= " + str(int(instance.t_var.value)) + "\n")
        f.write("CONSTRAINTS\n" + "".join(cons) + "\n")

        res = ["Selected epitopes\t" + str(len(result)) + ""]

        if int(instance.t_var.value) > 0:
            cov_anti = []
            for an in instance.Q:
                for e in result:
                    if e in instance.E_var[an].value:
                        cov_anti.append(an)
            cov_anti = set(cov_anti)
            res.append("Covered antigens\t" + str(len(cov_anti)) + " of " + str(len(instance.Q)) + "")
        cov_als = []
        res_set = set(result)
        locus = {}
        for a in instance.A:
            eps_of_all_i = list(instance.A_I[a])
            if res_set.intersection(set(eps_of_all_i)):
                cov_als.append(a)
                locus.setdefault(str(a).split("*")[0], set()).add(a)
        cov_als = set(cov_als)
        res.append("Covered alleles\t" + str(len(cov_als)) + " of " + str(len(instance.A)) + "")
        res.append("Locus coverage:\n")

        pop_cov = 1
        for k, g in locus.iteritems():
            locus = list(g)
            pop_cov *= (1.0 - sum(float(instance.p[a]) for a in locus)) ** 2
            covered = len(locus) / sum(1 for a in instance.A if a.split("*")[0] == k)
            res.append("\t%s\t%.2f" % (k, covered * 100))
        res.append("Population coverage:\t\t%.2f" % ((1.0 - pop_cov) * 100))
        f.write("RESULTS\n" + "\n".join(res) + "\n\n")

        is_antigen_cons = int(instance.t_var.value) > 0
        header = "Epitope\tConservation\tFraction of overall immunogenicity\tCovered alleles%s\n" % (
            "\tCovered antigens" if is_antigen_cons else "")

        rows = []
        overall_imm = sum(float(instance.i[e, a]) * float(instance.p[a]) for e in result for a in instance.A)
        for e in result:
            row = str(e) + "\t"
            if float(instance.t_c.value) > 0:
                row += str(float(instance.c[e]) * 100) + "\t"
            else:
                row += "100%\t"
            row += "%0.2f\t" % (sum(float(instance.i[e, a]) * float(instance.p[a]) for a in instance.A) / overall_imm)
            row += "%s" % " ".join(str(a) for a in instance.A if e in instance.A_I[a])
            # row += "%s"%"\n".join( str(float(instance.i[e,a])*float(instance.p[a]))  for a in instance.A if e in instance.A_I[a] )
            if is_antigen_cons:
                row += "\t%s" % " ".join(str(q) for q in instance.Q if e in instance.E_var[q])
            # row+="\n"
            rows.append(row)


        f.write(header + "\n".join(rows) + "\n\n")
        f.write("TARGET POPULATION/INDIVIDUAL\n")
        f.write("\tAllele\tProbability\n")
        for a in instance.A:
            f.write("\t%s\t%.5f\n" % (str(a), float(instance.p[a])))


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

    parser.add("cons_allele",
               required=False,
               type=float,
               default=None,
               help="Activates allele coverage constraint with specified threshold",
               short_name="cal")
    parser.add("cons_antigen",
               required=False,
               type=float,
               default=None,
               help="Activates antigen coverage constraint with specified threshold",
               short_name="can")
    parser.add("cons_conservation",
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
    thresh = {a.name: float(args["threshold"]) for a in epitopePrediciton.columns}
    opti = OptiTope(epitopePrediciton, threshold=thresh, k=int(args["k"]), solver="cplex", verbosity=0)

    # set constraints
    if args["cons_allele"] is not None:
        #print "allele constraint enforced"
        opti.activate_allele_coverage_const(float(args["cons_allele"]) / 100)

    if args["cons_antigen"] is not None:
        opti.activate_antigen_coverage_const(float(args["cons_antigen"]) / 100)

    if args["cons_conservation"] is not None:

        if args["conservation"] is not None:
            conservation = {}
            with open(args["conservation"], "r") as f:
                for l in f:
                    seq, cons = l.replace(",", " ").replace(";", " ").split()
                    conservation[seq.strip().upper()] = float(cons.strip())
            opti.activate_epitope_conservation_const(float(args["cons_conservation"])/100.0, conservation=conservation)
        else:
            opti.activate_epitope_conservation_const(float(args["cons_conservation"])/100.0)
    try:
        result = opti.solve(options={"threads":1})
        print "Result: ", result
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