#script for predicting snp's effect on TF binding using atSNP

#construct a specific directory
datadir=testdata

#1 modify names of TFs atSNP can predict as their standard HGNC SYMBOL

#JASPAR 2022
wget -c https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt --no-check-certificate
grep '>' JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt > jaspar2022.name.txt
sed -i 's/>//g' jaspar2022.name.txt
sed -i '1i\motif_id\ttf' jaspar2022.name.txt

#HOCOMOCO v11
wget -c https://hocomoco11.autosome.org/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_HUMAN_mono_jaspar_format.txt --no-check-certificate
wget -c https://hocomoco11.autosome.org/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv --no-check-certificate

Rscript TFrename.R


#2 select some TFs to predict (TFs with >= 20 pbsnps for predicting)
mkdir $datadir
Rscript testdata.predictedTF.R #note: you need to change the code for your data


#3 model's input snp: snpid chr snp(position) a1(ref allele) a2(alt allele)
Rscript --vanilla snp.R -f ../snpdata/testdata/testsnppos.tsv -o $datadir/testsnppos.tsv
#input for the script: -f:a file of snps,a snp per line,eg:chr10_114258723_G_A


#4 model's prediction

#4.1 JASPAR 2022
mkdir $datadir/pwm_jaspar2022
mkdir $datadir/pwm_jaspar2022/results
nohup Rscript --vanilla atSNP.predict.R -f $datadir/testsnppos.tsv -m $datadir/evaldata_interatSNP.jaspar587tf.csv --pwmdb jaspar -g BSgenome.Hsapiens.UCSC.hg19 -n 20 -o $datadir/pwm_jaspar2022/results/atSNP.jaspar2022.txt &
#input for the script: -f:output of the script 'snp.R', a file including 5 columns:snpid chr snp(position) a1(ref allele) a2(alt allele)
#-m: TF model you need to run,a TF model per line. The 'motif_name' column must be provided,eg:MA0069.1
#output for the script: see https://www.bioconductor.org/packages/release/bioc/vignettes/atSNP/inst/doc/atsnp-vignette.html

#4.2 HOCOMOCO v11
mkdir $datadir/pwm_hocomocov11
mkdir $datadir/pwm_hocomocov11/results
nohup Rscript --vanilla atSNP.predict.R -f $datadir/testsnppos.tsv -m $datadir/evaldata_interatSNP.hocomoco401tf.csv --pwmdb hocomoco -g BSgenome.Hsapiens.UCSC.hg19 -n 20 -o $datadir/pwm_hocomocov11/results/atSNP.hocomocov11.txt
#-m: TF model you need to run,a TF model per line. The 'motif_name' column must be provided,eg:ALX1_HUMAN.H11MO.0.B


#5 merge experimental and predictive difference value(2 alleles off snp) of TF binding
Rscript testdata.merge.R #note: you need to change the code for your data

#6 calculate AUROC,AUPRC of TFs

#6.1 JASPAR 2022
Rscript --vanilla auroc_auprc.R -e ../snpdata/testdata/evaldata -f $datadir/atSNP.jaspar2022.merged.expe.pred.results.txt -m $datadir/evaldata_interatSNP.jaspar587tf.csv -d log_lik_ratio -o $datadir/atSNP.jaspar2022.tf.roc.prc.txt
mv besttfmodel.roc.prc.txt atSNP.jaspar2022.besttfmodel.roc.prc.txt
#input for the script: -e: prefix of positive set and negative set,eg: evaldata_positive_data.txt,evaldata_negative_data.txt .The 'snp' and 'TF_SYMBOL' columns must be provided.
#-f: a file of merged experimental and predictive difference value(2 alleles of snp) of TF binding, the 'snp','TF_SYMBOL','model_name',predictive difference value of TF binding colums must be provided.
#-d: predictive difference value(2 alleles of snp) of TF binding

#6.2 HOCOMOCO v11
Rscript --vanilla auroc_auprc.R -e ../snpdata/testdata/evaldata -f $datadir/atSNP.hocomocov11.merged.expe.pred.results.txt -m $datadir/evaldata_interatSNP.hocomoco401tf.csv -d log_lik_ratio -o $datadir/atSNP.hocomocov11.tf.roc.prc.txt

# 7 regardless of pwm source,select TF's max AUROC and best models' merged experimental and predictive difference value(2 alleles off snp) of TF binding per TF,the data is used for computing correlation
Rscript --vanilla bestmotif.R --f1 $datadir/atSNP.jaspar2022.tf.roc.prc.txt --f2 $datadir/atSNP.hocomocov11.tf.roc.prc.txt --f3 $datadir/atSNP.jaspar2022.merged.expe.pred.results.txt --f4 $datadir/atSNP.hocomocov11.merged.expe.pred.results.txt --o1 $datadir/atSNP.alltf.bestmodel.roc.prc.txt --o2 $datadir/atSNP.alltf.merged.expe.pred.results.txt

