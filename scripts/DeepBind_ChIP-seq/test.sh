#script for predicting snp's effect on TF binding using DeepBind_ChIP-seq

#construct a specific directory
datadir=testdata

#1 DeepBind from kipoi
kipoi env create DeepBind

#2 modify names of TFs DeepBind_ChIP-seq can predict as their standard HGNC SYMBOL
awk -F '\t' '{print $1}' ~/.kipoi/models/DeepBind/models.tsv|grep 'Homo_sapiens/TF/'|grep 'ChIP-seq' > hg.chipseq.137tf.models.txt
Rscript TFrename.R

#3 select some TFs to predict (TFs with >= 20 pbsnps for predicting)
mkdir $datadir
Rscript testdata.predictedTF.R #note: you need to change the code for your data

#4 model's input sequence: 101-bp
python generate_allelic_seqs.py -f ../../genome/hg19.fa -s ../../snpdata/testdata/testsnppos.tsv -o $datadir/testdata 2>$datadir/testsnp.log
#input for the script: -f:reference genome; -s:a file of snps,a snp per line,eg:chr10_114258723_G_A

#5 model's prediction
mkdir $datadir/results
conda activate kipoi-DeepBind
nohup python deepbind.predict.py -f $datadir/testdata -m $datadir/evaldata_interchipseq136tf.csv -o $datadir/results/ &
#input for the script: -f:output of the script 'generate_allelic_seqs.py',prefix of the files including the sequences with ref or alt allele;
#-m: TF model you need to run,a TF model per line. The 'model_name' column must be provided,eg:DeepBind/Homo_sapiens/TF/D00317.009_ChIP-seq_CEBPB
#output for the script: files including 5 columns: snp, ref_pred, alt_pred, delta_alt_ref (alt_pred-ref_pred), TF(maybe non standard HGNC SYMBOL)

#6 merge experimental and predictive difference value(2 alleles off snp) of TF binding
Rscript testdata.merge.R #note: you need to change the code for your data

#7 calculate AUROC,AUPRC of TFs
Rscript --vanilla auroc_auprc.R -e ../../snpdata/testdata/evaldata -f $datadir/DeepBind_ChIP-seq.merged.expe.pred.results.txt -m $datadir/evaldata_interchipseq136tf.csv -d delta_alt_ref -o $datadir/DeepBind_ChIP-seq.tf.roc.prc.txt
#input for the script: -e: prefix of positive set and negative set,eg: evaldata_positive_data.txt,evaldata_negative_data.txt .The 'snp' and 'TF_SYMBOL' columns must be provided.
#-f: a file of merged experimental and predictive difference value(2 alleles of snp) of TF binding, the 'snp','TF_SYMBOL','model_name',predictive difference value of TF binding colums must be provided.
#-d: predictive difference value(2 alleles of snp) of TF binding


