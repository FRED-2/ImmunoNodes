#!/usr/bin/env python
"""
Command line tool for epitope selection

usage: epitopeselection.py [-h] -i INPUT -a ALLELES [-k K] [-t THRESHOLD] -o
                           OUTPUT [-s SOLVER] [-c_al CONS_ALLELE]
                           [-c_a CONS_ANTIGEN] [-c_c CONS_CONSERVATION]
                           [-c CONSERVATION]

Epitope Selection for vaccine design.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Peptide with immunogenicity file (from
                        epitopeprediction)
  -a ALLELES, --alleles ALLELES
                        Allele file with frequencies (one allele and frequency
                        per line)
  -k K, --k K           Specifies the number of epitopes to select
  -t THRESHOLD, --threshold THRESHOLD
                        Specifies the binding threshold for all alleles
  -o OUTPUT, --output OUTPUT
                        Specifies the output path. Results will be written to
                        CSV
  -s SOLVER, --solver SOLVER
                        Specifies the ILP solver
  -c_al CONS_ALLELE, --cons_allele CONS_ALLELE
                        Activates allele coverage constraint with specified
                        threshold
  -c_a CONS_ANTIGEN, --cons_antigen CONS_ANTIGEN
                        Activates antigen coverage constraint with specified
                        threshold
  -c_c CONS_CONSERVATION, --cons_conservation CONS_CONSERVATION
                        Activates conservation constraint with specified
                        threshold
  -c CONSERVATION, --conservation CONSERVATION
                        Specifies a Conservation file. First column is the
                        peptide seq second column the conservation.

"""
import sys
import pandas
import collections
import argparse

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
            protPos = {Protein("", gene_id=p, transcript_id=p): [0] for p in str(row["Protein ID"]).split(",")}
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
        f.write("#Prediction method: " + pred_method + "\n#\n")

        cons = ["#Maximum number of epitopes to select = " + str(int(instance.k.value)) + "\n"]
        if float(instance.t_c.value) > 0:
            cons.append("#Epitope conservation >= " + str(float(instance.t_c.value) * 100) + "%\n")

        if float(instance.t_allele.value) > 0:
            cons.append("#Covered alleles >= " + str(int(instance.t_allele.value)) + "\n")

        if float(instance.t_var.value) > 0:
            cons.append("#Covered antigens >= " + str(int(instance.t_var.value)) + "\n")
        f.write("#CONSTRAINTS\n" + "".join(cons) + "#\n")

        res = ["#Selected epitopes\t" + str(len(result)) + ""]

        if int(instance.t_var.value) > 0:
            cov_anti = []
            for an in instance.Q:
                for e in result:
                    if e in instance.E_var[an].value:
                        cov_anti.append(an)
            cov_anti = set(cov_anti)
            res.append("#Covered antigens\t" + str(len(cov_anti)) + " of " + str(len(instance.Q)) + "")
        cov_als = []
        res_set = set(result)
        locus = {}
        for a in instance.A:
            eps_of_all_i = list(instance.A_I[a])
            if res_set.intersection(set(eps_of_all_i)):
                cov_als.append(a)
                locus.setdefault(str(a).split("*")[0], set()).add(a)
        cov_als = set(cov_als)
        res.append("#Covered alleles\t" + str(len(cov_als)) + " of " + str(len(instance.A)) + "")
        res.append("#Locus coverage:")

        pop_cov = 1
        for k, g in locus.iteritems():
            locus = list(g)
            pop_cov *= (1.0 - sum(float(instance.p[a]) for a in locus)) ** 2
            covered = len(locus) / sum(1 for a in instance.A if a.split("*")[0] == k)
            res.append("#\t%s\t%.2f" % (k, covered * 100))
        res.append("#Population coverage:\t\t%.2f" % ((1.0 - pop_cov) * 100))
        f.write("#RESULTS\n" + "\n".join(res) + "\n")

        is_antigen_cons = int(instance.t_var.value) > 0
        header = "Epitope\tConservation\tFraction of overall immunogenicity\tCovered alleles%s\n" % (
            "\tCovered antigens" if is_antigen_cons else "")

        rows = []
        overall_imm = sum(float(instance.i[e, a]) * float(instance.p[a]) for e in result for a in instance.A)
        for e in result:
            row = str(e) + "\t"
            if float(instance.t_c.value) > 0:
                row += str(float(instance.c[e].value) * 100) + "\t"
            else:
                row += "100%\t"
            row += "%0.2f\t" % (sum(float(instance.i[e, a]) * float(instance.p[a]) for a in instance.A) / overall_imm)
            row += "%s" % " ".join(str(a) for a in instance.A if e in instance.A_I[a])
            if is_antigen_cons:
                row += "\t%s" % " ".join(str(q) for q in instance.Q if e in instance.E_var[q])
            rows.append(row)

        f.write(header + "\n".join(rows) + "\n\n")


def main():
    '''
        some input stuff
    '''
    parser = argparse.ArgumentParser(
                      description="Epitope Selection for vaccine design.",

    )
    parser.add_argument("-i","--input",
               required=True,
               type=str,
               help="Peptide with immunogenicity file (from epitopeprediction)",
    )
    parser.add_argument("-a","--alleles",
               required=True,
               type=str,
               help="Allele file with frequencies (one allele and frequency per line)",
    )

    parser.add_argument("-k","--k",
               required=False,
               type=int,
               default=10,
               help="Specifies the number of epitopes to select",
    )
    parser.add_argument("-t", "--threshold",
               type=float,
               default=0.,
               help="Specifies the binding threshold for all alleles",
    )
    parser.add_argument("-o", "--output",
               required=True,
               type=str,
               help="Specifies the output path. Results will be written to CSV",
    )

    parser.add_argument("-s","--solver",
               type=str,
               default="cbc",
               help="Specifies the ILP solver")

    parser.add_argument("-c_al", "--cons_allele",
               required=False,
               type=float,
               default=0.0,
               help="Activates allele coverage constraint with specified threshold",
    )

    parser.add_argument("-c_a", "--cons_antigen",
               required=False,
               type=float,
               default=0.0,
               help="Activates antigen coverage constraint with specified threshold",
    )

    c_c = parser.add_argument("-c_c", "--cons_conservation",
               required=False,
               type=float,
               help="Activates conservation constraint with specified threshold",
    )

    parser.add_argument("-c", "--conservation",
               required=False,
               type=str,
               help="Specifies a Conservation file. First column is the peptide seq second column the conservation.",
               )

    args = parser.parse_args()

    epitopePrediciton, method = generate_epitope_result(args.input, args.alleles)
    thresh = {a.name: float(args.threshold) for a in epitopePrediciton.columns}
    opti = OptiTope(epitopePrediciton, threshold=thresh, k=int(args.k), solver=args.solver, verbosity=0)

    # set constraints
    if args.cons_allele > 0:
        #print "allele constraint enforced"
        opti.activate_allele_coverage_const(float(args.cons_allele) / 100)

    if args.cons_antigen > 0:
        opti.activate_antigen_coverage_const(float(args.cons_antigen) / 100)

    if args.cons_conservation > 0:
        if args.conservation:
            conservation = {}
            with open(args.conservation, "r") as f:
                for l in f:
                    if l != "":
                        seq, cons = l.replace(",", " ").replace(";", " ").split()
                        conservation[seq.strip().upper()] = float(cons.strip())
            opti.activate_epitope_conservation_const(float(args.cons_conservation)/100.0, conservation=conservation)
        else:
            opti.activate_epitope_conservation_const(float(args.cons_conservation)/100.0)
    try:
        result = opti.solve(options={"threads": 1})

        to_csv(args.output, result, opti.instance, method)
        return 0
    except ValueError as e:
        sys.stderr.write("Could not optimally solve the problem. Please modify your constraints.\n"+str(e))
        return -1
    except Exception as e:
        sys.stderr.write(str(e))
        return -1


if __name__ == "__main__":
    sys.exit(main())