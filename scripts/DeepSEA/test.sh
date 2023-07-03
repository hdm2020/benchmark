#script for predicting snp's effect on TF binding using DeepSEA

#construct a specific directory
datadir=testdata

#1 DeepSEA from kipoi
kipoi env create DeepSEA/predict

#2 modify names of TFs DeepSEA can predict as their standard HGNC SYMBOL
cp ~/.kipoi/models/DeepSEA/predictor.names ./
grep 'DNase\|ac|\|me1|\|me2|\|me3|\|H2AZ' predictor.names > nonTF.list
Rscript TFrename.R

#3 select some TFs to predict (TFs with >= 20 pbsnps for predicting)
mkdir $datadir
Rscript testdata.predictedTF.R #note: you need to change the code for your data

#4 model's input sequence: 1000-bp
python generate_allelic_seqs.py -f ../../genome/hg19.fa -s ../../snpdata/testdata/testsnppos.tsv -o $datadir/testdata 2>$datadir/testsnp.log
#input for the script: -f:reference genome; -s:a file of snps,a snp per line,eg:chr10_114258723_G_A

#5 model's prediction
mkdir $datadir/results
conda activate kipoi-DeepSEA__predict
nohup python deepsea.predict.py -f $datadir/testdata -o $datadir/results/DeepSEA.txt &
#input for the script: -f:output of the script 'generate_allelic_seqs.py',prefix of the files including the sequences with ref or alt allele;
#output for the script: a file including 920 columns(no header):snp & 919 features(difference value:alt-ref)

#6 merge experimental and predictive difference value(2 alleles of snp) of TF binding
Rscript testdata.merge.R #note: you need to change the code for your data

#7 calculate AUROC,AUPRC of TFs
Rscript --vanilla auroc_auprc.R -e ../../snpdata/testdata/evaldata -f $datadir/DeepSEA.merged.expe.pred.results.txt -m $datadir/evaldata_interdeepsea689tf.csv -d delta_alt_ref -o $datadir/DeepSEA.tf.roc.prc.txt
mv besttfmodel.roc.prc.txt $datadir/DeepSEA.alltf.bestmodel.roc.prc.txt
#input for the script: -e: prefix of positive set and negative set,eg: evaldata_positive_data.txt,evaldata_negative_data.txt .The 'snp' and 'TF_SYMBOL' columns must be provided.
#-f: a file of merged experimental and predictive difference value(2 alleles of snp) of TF binding, the 'snp','TF_SYMBOL','model_name',predictive difference value of TF binding colums must be provided.
#-d: predictive difference value(2 alleles of snp) of TF binding

#7 best models' merged experimental and predictive difference value(2 alleles off snp) of TF binding per TF,the data is used for computing correlation
Rscript --vanilla bestmodel.R --f1 $datadir/DeepSEA.alltf.bestmodel.roc.prc.txt --f2 $datadir/DeepSEA.merged.expe.pred.results.txt --out $datadir/DeepSEA.alltf.merged.expe.pred.results.txt

