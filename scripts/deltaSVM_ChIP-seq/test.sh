#script for predicting snp's effect on TF binding using deltaSVM_ChIP-seq

#construct a specific directory
datadir=testdata

#1 download deltaSVM
wget -c http://www.beerlab.org/deltasvm_models/downloads/deltasvm_models_e3tf.tar.gz --no-check-certificate
mkdir e3tfmodels
tar -zxvf deltasvm_models_e3tf.tar.gz -C e3tfmodels/
rm deltasvm_models_e3tf.tar.gz
#script for predicting SNP's effect
wget -c http://www.beerlab.org/deltasvm_models/downloads/score_snp_seq.py --no-check-certificate
#modify this script to predict SNP's effect:deltaSVM_ChIP-seq.predict.py

#2 modify names of TFs deltaSVM_ChIP-seq can predict as their standard HGNC SYMBOL
wget -c http://www.beerlab.org/deltasvm_models/downloads/sample_list.txt --no-check-certificate
grep TF_E3 sample_list.txt >tfe3.699.txt
ls e3tfmodels/ > 699.models.txt

Rscript TFrename.R

#3 select some TFs to predict (TFs with >= 20 pbsnps for predicting)
mkdir $datadir
Rscript testdata.predictedTF.R #note: you need to change the code for your data

#4 model's input sequence: 21-bp
python generate_allelic_seqs.py -f ../../genome/hg19.fa -s ../../snpdata/testdata/testsnppos.tsv -o $datadir/testdata 2>$datadir/testsnp.log
#input for the script: -f:reference genome; -s:a file of snps,a snp per line,eg:chr10_114258723_G_A

#4 model's prediction
mkdir $datadir/results
for weight in `cat testdata/predicted.tf.model.txt`
do
python deltaSVM_ChIP-seq.predict.py $datadir/testdata.ref.fa $datadir/testdata.alt.fa e3tfmodels/$weight $datadir/results/
done
#input for the script: output of the script 'generate_allelic_seqs.py',the sequences including ref or alt allele;
#weight file: a weight file per line, eg:TF_E3_104_hg38_300_top10k_vs_neg1x_avg_weights.out
#output for the script: files including 2 columns(no header): snp & difference value(alt-ref)

#5 merge experimental and predictive difference value(2 alleles off snp) of TF binding
Rscript testdata.merge.R #note: you need to change the code for your data

#6 calculate AUROC,AUPRC of TFs
Rscript --vanilla auroc_auprc.R -e ../../snpdata/testdata/evaldata -f $datadir/deltaSVM_ChIP-seq.merged.expe.pred.results.txt -m $datadir/evaldata_inter699e3tf.csv -d deltaSVM -o $datadir/deltaSVM_ChIP-seq.tf.roc.prc.txt
mv besttfmodel.roc.prc.txt $datadir/deltaSVM_ChIP-seq.alltf.bestmodel.roc.prc.txt
#input for the script: -e: prefix of positive set and negative set,eg: evaldata_positive_data.txt,evaldata_negative_data.txt .The 'snp' and 'TF_SYMBOL' columns must be provided.
#-f: a file of merged experimental and predictive difference value(2 alleles of snp) of TF binding, the 'snp','TF_SYMBOL','model_name',predictive difference value of TF binding colums must be provided.
#-d: predictive difference value(2 alleles of snp) of TF binding

#7 best models' merged experimental and predictive difference value(2 alleles off snp) of TF binding per TF,the data is used for computing correlation
Rscript --vanilla bestmodel.R --f1 $datadir/deltaSVM_ChIP-seq.alltf.bestmodel.roc.prc.txt --f2 $datadir/deltaSVM_ChIP-seq.merged.expe.pred.results.txt --out $datadir/deltaSVM_ChIP-seq.alltf.merged.expe.pred.results.txt


