#script for predicting snp's effect on TF binding using DeepSEA

import sys
import os
import numpy as np
import pandas as pd
import kipoi
import argparse

#convert sequences into one-hot encoded format
def seq_to_hot(seq):
#lowercase to uppercase
    seq=seq.replace('a','A')
    seq=seq.replace('c','C')
    seq=seq.replace('g','G')
    seq=seq.replace('t','T')
    seq=seq.replace('n','N')
#encode A:1, other bases:0
    Aseq=seq
    Aseq=Aseq.replace('A','1')
    Aseq=Aseq.replace('C','0')
    Aseq=Aseq.replace('G','0')
    Aseq=Aseq.replace('T','0')
    Aseq=Aseq.replace('N','0')
    Aseq=np.asarray(list(Aseq),dtype='float32')
#encode C:1, other bases:0
    Cseq=seq
    Cseq=Cseq.replace('A','0')
    Cseq=Cseq.replace('C','1')
    Cseq=Cseq.replace('G','0')
    Cseq=Cseq.replace('T','0')
    Cseq=Cseq.replace('N','0')
    Cseq=np.asarray(list(Cseq),dtype='float32')
#encode G:1, other bases:0
    Gseq=seq
    Gseq=Gseq.replace('A','0')
    Gseq=Gseq.replace('C','0')
    Gseq=Gseq.replace('G','1')
    Gseq=Gseq.replace('T','0')
    Gseq=Gseq.replace('N','0')
    Gseq=np.asarray(list(Gseq),dtype='float32')
#encode T:1, other bases:0
    Tseq=seq
    Tseq=Tseq.replace('A','0')
    Tseq=Tseq.replace('C','0')
    Tseq=Tseq.replace('G','0')
    Tseq=Tseq.replace('T','1')
    Tseq=Tseq.replace('N','0')
    Tseq=np.asarray(list(Tseq),dtype='float32')
#one-hot encoded sequences
    hot=np.vstack((Aseq,Cseq,Gseq,Tseq))
#convert format: (4，101) -- (1,4,1,1000)
    x=np.array([hot])
    y=np.moveaxis(x,0,1)
    hot=np.array([y])  #4 dimension
    return hot

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

#prediction
def deepsea_predict(ref_file,alt_file):
    ref_seq_dict,alt_seq_dict=seq_to_dict(ref_file=ref_file,alt_file=alt_file)
    model = kipoi.get_model('DeepSEA/predict')
    df=[]
    df=pd.DataFrame(df)
    for pos_name in ref_seq_dict.keys():
        ref_seq=ref_seq_dict[pos_name]
        alt_seq=alt_seq_dict[pos_name]
        ref_seq_onehot=seq_to_hot(ref_seq)
        alt_seq_onehot=seq_to_hot(alt_seq)
        ref_pred=model.predict_on_batch(ref_seq_onehot)
        alt_pred=model.predict_on_batch(alt_seq_onehot)
        diff=alt_pred-ref_pred
        diff=pd.DataFrame(diff)
        diff.insert(loc=0, value=pos_name,column='snp')
        df=pd.concat([df,diff])
    return df

#
def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-f", "--fasta", dest="fasta", required=True, help="oligo fasta")
    parser.add_argument("-o", "--out", dest="out", required=True, help="output file/directory (bindind difference)")
    args = parser.parse_args()
    
    #sequences
    reffa = args.fasta+".ref.fa"
    altfa = args.fasta+".alt.fa"
    
    #output
    outfile = args.out

    #predict
    batch_result = deepsea_predict(ref_file=reffa,alt_file=altfa)
    
    #save results
    batch_result.to_csv(outfile,index=False,header=False,sep='\t')


if __name__ == "__main__":
    sys.exit(main())

