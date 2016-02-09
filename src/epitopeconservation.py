#!/usr/bin/env python
"""
Commandline tool for epitope prediction

usage: epitopeconservation.py [-h] -i INPUT [-l LENGTH] -cons
                              OUTPUT_CONSERVATION [-f OUTPUT_FASTA]

Epitope conservation computation method.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Special MSA file
  -l LENGTH, --length LENGTH
                        Specifies the length of the peptides (default=9).
  -cons OUTPUT_CONSERVATION, --output_conservation OUTPUT_CONSERVATION
                        specifies output file
  -f OUTPUT_FASTA, --output_fasta OUTPUT_FASTA
                        specifies output file

"""
from __future__ import division
import numpy
import sys
import argparse
from Fred2.Core.Protein import Protein
from Fred2.Core.Peptide import Peptide


def isValidAASequence(epitope):
    epitope = epitope.strip()
    if epitope == '':
        return 1

    if not epitope.isalpha():
        return 0

    # check for invalid letters
    for aa in ['B', 'J', 'O', 'U', 'X', 'Z']:
        if aa in epitope:
            return 0

    return 1


def isValidMSASequence(seq):
    seq = seq.upper()

    return isValidAASequence(seq.replace('-', ''))


def extractEpitopesAndConservationFromConsensus(consensus_info, epitope_length):
    # print "consensus_info", consensus_info
    error = ''
    epitope_data = []
    epitope_antigen_data = {}
    antigens = consensus_info.keys()
    conservation = {}
    epitopes = {}

    # print "type of consensus_info ",type(consensus_info)
    for ci in consensus_info.values():
        consensus = str(ci[0][0])
        frequencies = ci[0][1]

        # consensus = consensus.upper()
        # if isValidAASequence(consensus):

        for i in xrange(len(consensus) - epitope_length + 1):
            epitope = Peptide(consensus[i:i + epitope_length])
            # print "peptide", epitope
            co = numpy.product(frequencies[i:i + epitope_length])

            if not conservation.has_key(epitope):
                # print "test type prot in epitopes", epitope.proteins
                conservation[epitope] = co
                epitopes[epitope] = epitope
            else:
                epitope = epitopes[epitope]
                # print "test type prot in epitopes", epitope.proteins
                if conservation[epitope] < co:
                    conservation[epitope] = co
    return (error, conservation)


def determineConsensusFromMSA(sequences):
    consensus_sequence = ''
    frequencies = numpy.array([])
    error = ''
    aa_tracker = []

    # print sequences
    number_of_sequences = len(sequences)
    if number_of_sequences > 0:
        seq_len = len(sequences[0])
        frequencies = numpy.zeros(seq_len)

        for i in xrange(seq_len):
            aa_tracker.append({})

        for seq in sequences:
            if len(seq) == seq_len:

                for p in xrange(seq_len):
                    c = seq[p]
                    if aa_tracker[p].has_key(c):
                        aa_tracker[p][c] += 1
                    else:
                        aa_tracker[p][c] = 1

            else:  # sequences of different lengths
                error = "MSA sequences are of different lengths."
                break

        if error == '':
            for i in range(len(aa_tracker)):
                max = 0
                max_aa = ''
                for aa, count in aa_tracker[i].items():
                    if count > max:
                        max = count
                        max_aa = aa

                if max_aa != '-':  # don't include gaps
                    consensus_sequence += max_aa
                    frequencies[i] = 1. * max / number_of_sequences

    return consensus_sequence, frequencies


def extractEpitopeInformationFromMSA(input, epitope_length):
    error = ''
    NONE = 0
    ANTIGEN_HEADER = 1
    SEQUENCE_HEADER = 2
    SEQUENCE = 3

    last_read = NONE
    antigen = ''
    current_sequence = ''
    sequences = []

    consensus_info = {}
    consensuses = {}
    epitope_info = []
    antigens = []
    conservation = {}
    with open(input, "r") as f:
        for line in f:
            line = line.strip()
            if line == '' or line[0] == '#':
                continue
            # print 'in extract msa info',line
            if line.startswith('>'):
                if last_read == NONE or last_read == SEQUENCE:
                    if last_read == SEQUENCE:  # done with previous antigen

                        sequences.append(current_sequence)
                        current_sequence = ''

                    antigen = line[2:].strip()
                    last_read = SEQUENCE_HEADER

                else:  # wrong format
                    error = "Input format does not comply with required MSA format."
                    break

            elif line[0] == '>' and antigen != '':
                if last_read == ANTIGEN_HEADER or last_read == SEQUENCE:
                    if last_read == SEQUENCE:
                        sequences.append(current_sequence)
                        current_sequence = ''

                    last_read = SEQUENCE_HEADER

            else:  # sequence

                if not isValidMSASequence(line) and antigen != '':
                    error = "Invalid amino acid sequence given for antigen " + antigen + ". "
                    break

                if last_read == SEQUENCE_HEADER or last_read == SEQUENCE:
                    # current_sequence += line[:-1]
                    current_sequence += line
                    last_read = SEQUENCE
                else:
                    error = "Input format does not comply with required MSA format."
                    break

    if error == '':
        if last_read == SEQUENCE:
            sequences.append(current_sequence)
            consensus = determineConsensusFromMSA(sequences)
            if consensus[0] != '':
                consensuses[antigen] = consensus[0].upper()
                consensus_info.setdefault(antigen, []).append(
                    (Protein(consensus[0].upper(), antigen, antigen), consensus[1]))
            else:  # different sequence lengths
                error = "MSA sequences of antigen " + antigen + " are of different lengths."
        else:
            error = "Input format does not comply with required MSA format."

        if error == '':
            (error, conservation) = extractEpitopesAndConservationFromConsensus(consensus_info, epitope_length)
        # print "infos form extractEpitopesAndConservationFromConsensus", conservation

    return (error, conservation, consensuses)


def main():
    '''
        some input stuff
    '''
    parser = argparse.ArgumentParser(description="Epitope conservation computation method.")
    parser.add_argument("-i", "--input",
                        required=True,
                        help="Special MSA file"
                        )
    parser.add_argument("-l", "--length",
                        required=False,
                        type=int,
                        default=9,
                        help="Specifies the length of the peptides (default=9).")
    parser.add_argument("-cons", "--output_conservation",
                        required=True,
                        help="specifies output file")
    parser.add_argument("-f", "--output_fasta",
                        required=False,
                        default="",
                        help="specifies output file")
    args = parser.parse_args()

    # print "DEBUG", args

    info = extractEpitopeInformationFromMSA(args.input, args.length)

    if info[0] != "":
        raise Exception(info[0])

    with open(args.output_conservation, "w") as f:
        f.write("".join("%s\t%0.3f\n" % (str(epi), cons) for epi, cons in info[1].iteritems()))

    with open(args.output_fasta, "w") as f:
        f.write("\n".join(">%s\n%s" % (k, v) for k, v in info[2].iteritems()))


if __name__ == "__main__":
    main()
