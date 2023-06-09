#script for predicting snp's effect on TF binding using deltaSVM_ChIP-seq,modify from score_snp_seq.py
import os, sys
import numpy as np
import pandas as pd

#sequence information: 'dict': keys:peak_name(eg:chr10_100027987_C_G); values: sequences
def seq_to_dict(ref_file,alt_file):
    ref_seq_dict={}
    for line in open(ref_file):
        line=line.strip()
        if line[0] == '>':
            peak_name=line[1:]
        else:
            ref_seq_dict[peak_name] = line.upper()
    alt_seq_dict={}
    for line in open(alt_file):
        line=line.strip()
        if line[0] == '>':
            peak_name=line[1:]
        else:
            alt_seq_dict[peak_name] = line.upper()
    return ref_seq_dict,alt_seq_dict

def invert_char(str):
    str=str.upper()
    if str=="A":
        return "T"
    elif str=="T":
        return "A"
    elif str=="C":
        return "G"
    elif str=="G":
        return "C"
    else:
        return str
    

def invert(str):
    nstr=""
    for c in str:
        nstr=invert_char(c)+nstr
    return nstr

def main(argv = sys.argv):
    if len(argv) != 5:
        print("Usage: python", argv[0],"<seq1> <seq2> <weightfile> <output:directory>")
        exit(0)

    ref_seq_dict,alt_seq_dict=seq_to_dict(ref_file=argv[1],alt_file=argv[2])
    wtfile = argv[3]

    infile=open(wtfile,"r")
    type=0
    wt={}
    for line in infile:
        f=line.strip().split()
        k=len(f[0])
        wt[f[0]]=f[1]
        wt[invert(f[0])]=f[1];
    infile.close()

#    if len(seq1)<2*k-1 or len(seq2)<2*k-1 :
#        print("seqs should be at least {0}bp to accurately score with {1}-mer weights".format(2*k-1,k))
#        exit()
    df={}
    for pos_name in ref_seq_dict.keys():
        ref_seq=ref_seq_dict[pos_name]
        seq1 = ref_seq.upper()
        alt_seq=alt_seq_dict[pos_name]
        seq2 = alt_seq.upper()
        s1=0
        s2=0
        for i in range(0,len(seq1)-k+1):
            a=seq1[i:i+k]
            s1=s1+float(wt[a])
            b=seq2[i:i+k]
            s2=s2+float(wt[b])
        deltaSVM="{0:6.2f}".format(s2-s1)
        df[pos_name]=deltaSVM
    diff=pd.DataFrame.from_dict(df,orient='index',columns=['deltaSVM'])
    diff=diff.reset_index().rename(columns={'index':'snp'})
    outfile=argv[4]+argv[3].split('/')[-1]
    diff.to_csv(outfile,index=False,header=False,sep='\t')
    print(argv[3].split('/')[-1]+'finished')

if __name__ == '__main__': main()

