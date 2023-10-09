#script for predicting snp's effect on TF binding using Sei

#construct a specific directory
datadir=testdata

#1 Download Sei

#dependices
#pytorch
conda install pytorch==1.9.0 torchvision==0.10.0 torchaudio==0.9.0 cpuonly -c pytorch
#selen-sdk
git clone https://github.com/FunctionLab/selene.git
cd selene
python setup.py build_ext --inplace
python setup.py install
cd ../
#docopt
conda install docopt
#sei
git clone  https://github.com/FunctionLab/sei-framework.git
mv sei-framework/* ./
rm -rf sei-framework/
#download trained Sei model and resources (containing hg19 and hg38 FASTA files)
sh ./download_data.sh

#2 modify names of TFs Sei can predict as their standard HGNC SYMBOL
#'model/target.names' includes 3 types of chromatin feaures: histone modification,chromatin accessibility,TF binding

#position index of histone marks('histone_inds.npy',0-based)
#!!!you need to run codes below.
#import numpy as np
#import pandas as pd
#histone_inds = np.load('model/histone_inds.npy')
#df=pd.DataFrame(histone_inds)
#df.to_csv('histone_inds.txt',index=False,header=False,sep='\t')

#chromatin accessibility
grep 'DNase\|ATAC-seq' model/target.names > chrom_access.txt

Rscript TFrename.R


#3 select some TFs to predict (TFs with >= 20 pbsnps for predicting)
mkdir $datadir
Rscript testdata.predictedTF.R #note: you need to change the code for your data

#4 model's input:
Rscript --vanilla snp.R -f ../../snpdata/testdata/testsnppos.tsv -o $datadir/testsnppos.vcf
#input for the script: -f:a file of snps,a snp per line,eg:chr10_114258723_G_A


#5 model's prediction
mkdir $datadir/results
nohup sh 1_variant_effect_prediction.sh $datadir/testsnppos.vcf hg19 $datadir/results/ &
#Example usage:sh 1_variant_effect_prediction.sh <vcf> <hg> <output-dir> [--cuda]


#6 merge experimental and predictive difference value(2 alleles of snp) of TF binding
Rscript testdata.merge.R #note: you need to change the code for your data


#7 calculate AUROC,AUPRC of TFs
Rscript --vanilla auroc_auprc.R -e ../../snpdata/testdata/evaldata -f $datadir/Sei.merged.expe.pred.results.txt -m $datadir/evaldata_interSei9315tf.csv -d delta_alt_ref -o $datadir/Sei.tf.roc.prc.txt
mv besttfmodel.roc.prc.txt $datadir/Sei.alltf.bestmodel.roc.prc.txt
#input for the script: -e: prefix of positive set and negative set,eg: evaldata_positive_data.txt,evaldata_negative_data.txt .The 'snp' and 'TF_SYMBOL' columns must be provided.
#-f: a file of merged experimental and predictive difference value(2 alleles of snp) of TF binding, the 'snp','TF_SYMBOL','model_name',predictive difference value of TF binding colums must be provided.
#-d: predictive difference value(2 alleles of snp) of TF binding


#8 best models' merged experimental and predictive difference value(2 alleles off snp) of TF binding per TF,the data is used for computing correlation
Rscript --vanilla bestmodel.R --f1 $datadir/Sei.alltf.bestmodel.roc.prc.txt --f2 $datadir/Sei.merged.expe.pred.results.txt --out $datadir/Sei.alltf.merged.expe.pred.results.txt

