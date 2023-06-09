#script for predicting snp's effect on TF binding using deltaSVM_ChIP-seq


#construct a specific directory
datadir=testdata

#1 download deltaSVM

#1.1 94 high-confidence TFs
cd ../
wget -c https://github.com/ren-lab/deltaSVM/archive/refs/heads/master.zip --no-check-certificate
unzip master.zip 
mv deltaSVM-master deltaSVM_HT-SELEX

#download resources files(including reference genome and weights)
#see https://github.com/ren-lab/deltaSVM/README.md,The requested URL is unvalid now.
#get one directory: resources(files:hg19.fa  hg19.fa.fai  models.weights.txt  threhsolds.obs.tsv  threhsolds.pbs.tsv)

#1.2 all gkm-SVM models
#see https://github.com/ren-lab/deltaSVM/README.md,The requested URL is unvalid now.
#two directory: 533gkmsvm_models,deltasvm_weights


#2 modify names of TFs deltaSVM_HT-SELEX can predict as their standard HGNC SYMBOL
#94 high-confidence TFs
ls gkmsvm_models/ > 94tf.model.txt
#all 533 TFs 
ls 533gkmsvm_models/ > 533gkmsvm_models.txt

Rscript TFrename.R


#3 select some TFs to predict (TFs with >= 20 pbsnps for predicting)
mkdir $datadir
Rscript testdata.predictedTF.R #note: you need to change the code for your data


#4 model's input sequence: 40-bp
python scripts/generate_allelic_seqs.py -f resources/hg19.fa -s ../snpdata/snppos.tsv -o $datadir/testdata 2>$datadir/testdata.log
#input for the script: -f:reference genome; -s:a file of snps,a snp per line,eg:chr10_114258723_G_A


#5  model's prediction

#5.1 94 high-confidence TFs
mkdir $datadir/94tf_results
scripts/deltasvm_subset_multi $datadir/testdata.ref.fa $datadir/testdata.alt.fa resources/models.weights.txt $datadir/94tf_results/testdata.pbs.pred.tsv resources/threhsolds.pbs.tsv 2>$datadir/94tf_results/testdata.deltasvm.log
sort -k1,1 -k3,3 -o $datadir/94tf_results/testdata.pbs.pred.tsv $datadir/94tf_results/testdata.pbs.pred.tsv
sed -i '1i snp\tdeltaSVM\tTF\tpreferred_allele' $datadir/94tf_results/testdata.pbs.pred.tsv

#5.2 other 439 TFs
mkdir $datadir/439tf_results
mkdir $datadir/439tf_log
for model in `cat $datadir/predicted.tf.model.txt`
do
perl scripts/deltasvm.pl $datadir/testdata.ref.fa $datadir/testdata.alt.fa 2>$datadir/439tf_log/${model}.deltasvm.log deltasvm_weights/${model}.weights.txt $datadir/439tf_results/${model}.deltasvm_scores.tsv
done


#6 merge experimental and predictive difference value(2 alleles off snp) of TF binding
Rscript testdata.merge.R #note: you need to change the code for your data

#7 calculate AUROC,AUPRC of TFs

Rscript --vanilla auroc_auprc.R -e ../snpdata/testdata/evaldata -f $datadir/deltaSVM_HT-SELEX.merged.expe.pred.results.txt -m $datadir/evaldata_inter533tf.csv -d deltaSVM -o $datadir/deltaSVM_HT-SELEX.tf.roc.prc.txt
mv $datadir/besttfmodel.roc.prc.txt $datadir/deltaSVM_HT-SELEX.alltf.bestmodel.roc.prc.txt
#input for the script: -e: prefix of positive set and negative set,eg: evaldata_positive_data.txt,evaldata_negative_data.txt .The 'snp' and 'TF_SYMBOL' columns must be provided.
#-f: a file of merged experimental and predictive difference value(2 alleles of snp) of TF binding, the 'snp','TF_SYMBOL','model_name',predictive difference value of TF binding colums must be provided.
#-d: predictive difference value(2 alleles of snp) of TF binding

#8 best models' merged experimental and predictive difference value(2 alleles off snp) of TF binding per TF,the data is used for computing correlation
Rscript --vanilla bestmodel.R --f1 $datadir/deltaSVM_HT-SELEX.alltf.bestmodel.roc.prc.txt --f2 $datadir/deltaSVM_HT-SELEX.merged.expe.pred.results.txt --out $datadir/deltaSVM_HT-SELEX.alltf.merged.expe.pred.results.txt


