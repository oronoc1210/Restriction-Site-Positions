#! /usr/bin/env python
# Generate site positions in genome from given restriction enzyme
# Juicer 1.5
# Modified by Conor O'Donoghue on 3/7/19
#   Added multiple restriction site and ambiguity functionality
import sys
import os
import argparse
import gzip
import subprocess
import re

def parser_gen():
    parser = argparse.ArgumentParser(
        description="""Generates a list of restriction enzyme cuts in the genome.
        Can calculate cut sites from multiple restriction enzymes,
        and can account for ambiguity 
        (https://www.bioinformatics.org/sms/iupac.html).""",
        epilog="""Example:
        python generate_site_positions_modified.py 
        --cut_site_sequences GATC CCRYGG --genome hg19 --filepath test.fa""")
    parser.add_argument("--cut_site_sequences", 
                        dest='cut_site_sequences', nargs='+',
                        help="""Cut site sequences (e.g. GATC).
                             Separate multiple sequences with a single space.
                             Allows for ambiguity (e.g. CCRYGG).""")
    parser.add_argument("--genome",
                        dest='genome',
                        help="""Genome name. Used in outfile name.""")
    parser.add_argument("--fasta",
                        dest='fasta',
                        help="""Path to fasta file. Unzips if gzipped.""")
    return parser

def check_args(parser):
    args = parser.parse_args()
    if not all((args.cut_site_sequences, args.genome, args.fasta)):
        print("All three parameters must be entered! Exiting...")
        parser.print_help()
        sys.exit(0)
    if args.fasta and not os.path.exists(args.fasta):
        print("File %s does not exist! Exiting..." % args.fasta)
        sys.exit(1)
    if args.fasta and not any(x in args.fasta for x in ['.fa', '.fasta']):
        print("Input file %s must be fasta format! Exiting..." % args.fasta)
        sys.exit(1)
    return args

def main():
    args = check_args(parser_gen())
    results_list = []
    iupac = {'A': ['A'],
             'C': ['C'],
             'G': ['G'],
             'T': ['T'],
             'R': ['A', 'G'],
             'Y': ['C', 'T'],
             'S': ['G', 'C'],
             'W': ['A', 'T'],
             'K': ['G', 'T'],
             'M': ['A', 'C'],
             'B': ['C', 'G', 'T'],
             'D': ['A', 'G', 'T'],
             'H': ['A', 'C', 'T'],
             'V': ['A', 'C', 'G'],
             'N': ['A', 'C', 'T', 'G']
             }
    if '.gz' in args.fasta:
        f=gzip.open(args.fasta, 'rb')
    else:
        f=open(args.fasta, 'r')
    g=open(args.genome+'_'+'_'.join(args.cut_site_sequences)+'.txt', 'w+')
    chromosomecheck=f.readline()
    while (chromosomecheck):
        if (chromosomecheck.startswith('>')):
            firststr=re.findall('(?<=>)([A-Za-z0-9_\- ]*)', chromosomecheck)[0]
            results_list.append(firststr + ' ')
            framesize = len(max(args.cut_site_sequences, key=len))
            nextchar=f.read(framesize)
            frame = []
            for i in nextchar:
                frame.append(i)
            position=1
            while len(nextchar) != 0:
                for cut_sequence in args.cut_site_sequences:
                    if all(frame[i] in iupac[cut_sequence[i]] for i in range(len(cut_sequence))):
                        results_list.append(str(position) + ' ')
                        break
                frame.pop(0)
                nextchar=f.read(1)
                if nextchar=='\n':
                    nextchar=f.read(1)
                if nextchar=='>':
                    break
                frame.append(nextchar)
                position+=1
            position+=len(frame) - 1
            results_list.append(str(position) + '\n')
            chromosomecheck=nextchar+f.readline()
        else:
            chromosomecheck=f.readline()
    print(chromosomecheck)
    g.write(''.join(results_list) + '\n')
    f.close()
    g.close()

if __name__ == '__main__':
    main()
