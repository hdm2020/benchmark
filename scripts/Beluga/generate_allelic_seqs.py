#!/usr/bin/python

#########################################
# Author: Yunjiang Qiu <serein927@gmail.com>
# File: generate_allelic_seqs.py
# Create Date: 2018-02-20 10:22:38
#########################################
#use this script from Yunjiang Qiu to generate 1000bp sequences,and 500th is reference allele or alternative allele of the SNP.
import sys
import argparse
import pysam

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-f", "--fasta", dest="fasta", required=True, help="oligo fasta")
    parser.add_argument("-s", "--snp", dest="snp", required=True, help="snp information")
    parser.add_argument("-o", "--out", dest="out", required=True, help="output fasta")
    args = parser.parse_args()

    fastafile = pysam.FastaFile(args.fasta)

    f1 = open(args.out+".ref.fa", "w")
    f2 = open(args.out+".alt.fa", "w")

    with open (args.snp, "r") as f:
        for line in f:
            try:
                chrom, pos, ref, alt = line.rstrip().split("_")
            except ValueError:
                continue
            ref = ref.upper()
            alt = alt.upper()
            if len(ref) > 1 or len(alt) > 1:
                continue
            oligo_seq = fastafile.fetch(chrom, int(pos)-1000, int(pos)+1000).upper()
            if oligo_seq[999] != ref and oligo_seq[999] != alt:
                sys.stderr.write("{0} is not matched to the genome.\n".format(line.rstrip()))
                continue
            f1.write(">{0}_{1}_{2}_{3}\n".format(chrom, pos, ref, alt))
            f1.write("{0}\n".format((oligo_seq[0:999] + ref + oligo_seq[1000:])))
            f2.write(">{0}_{1}_{2}_{3}\n".format(chrom, pos, ref, alt))
            f2.write("{0}\n".format((oligo_seq[0:999] + alt + oligo_seq[1000:])))

    f1.close()
    f2.close()

if __name__ == "__main__":
    sys.exit(main())

