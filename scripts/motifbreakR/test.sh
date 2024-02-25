#script for predicting snp's effect on TF binding using motifbreakR

#construct a specific directory
datadir=testdata

#1 modify names of TFs motifbrekR can predict as their standard HGNC SYMBOL
#JASPAR 2022 & HOCOMOCO v11

wget -c https://hocomoco11.autosome.org/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv --no-check-certificate

Rscript TFrename.R


#2 select some TFs to predict (TFs with >= 20 pbsnps for predicting)
mkdir $datadir
Rscript testdata.predictedTF.R #note: you need to change the code for your data


#3 model's input snp: bed format--chromosome, start, end, name, score and strand
Rscript --vanilla snp.R -f ../../snpdata/testdata/testsnppos.tsv -o $datadir/testsnppos.bed
#input for the script: -f:a file of snps,a snp per line,eg:chr10_114258723_G_A


#4 model's prediction

#JASPAR 2022
mkdir $datadir/pwm_jaspar2022
mkdir $datadir/pwm_jaspar2022/results
nohup Rscript --vanilla motifbreakR.predict.R -f $datadir/testsnppos.bed -m $datadir/evaldata_intermotifbreakR.jaspar587tf.csv --pwmdb jaspar -g BSgenome.Hsapiens.UCSC.hg19 --method default -o $datadir/pwm_jaspar2022/results/motifbreakR.jaspar2022.txt &
#input for the script: -f:output of the script 'snp.R', a file including 6 columns: chr start(position_start) snp(position_end) snpid score strand
#-m: TF model you need to run,a TF model per line. The 'motif_name' column must be provided,eg:Hsapiens-jaspar2022-FOXF2-MA0030.1
#output for the script: see https://bioconductor.org/packages/release/bioc/vignettes/motifbreakR/inst/doc/motifbreakR-vignette.html

#HOCOMOCO v11
mkdir $datadir/pwm_hocomocov11
mkdir $datadir/pwm_hocomocov11/results
nohup Rscript --vanilla motifbreakR.predict.R -f $datadir/testsnppos.bed -m $datadir/evaldata_intermotifbreakR.hocomoco400tf.csv --pwmdb hocomoco -g BSgenome.Hsapiens.UCSC.hg19 --method default -o $datadir/pwm_hocomocov11/results/motifbreakR.hocomocov11.txt &
#-m: TF model you need to run,a TF model per line. The 'motif_name' column must be provided,eg:AHR_HUMAN.H11MO.0.B

#5 merge experimental and predictive difference value(2 alleles off snp) of TF binding
Rscript testdata.merge.R #note: you need to change the code for your data

#6 calculate AUROC,AUPRC of TFs

#6.1 JASPAR 2022
Rscript --vanilla auroc_auprc.R -e ../../snpdata/testdata/evaldata -f $datadir/motifbreakR.jaspar2022.merged.expe.pred.results.txt -m $datadir/evaldata_intermotifbreakR.jaspar587tf.csv -d alleleDiff -o $datadir/motifbreakR.jaspar2022.tf.roc.prc.txt
mv besttfmodel.roc.prc.txt $datadir/motifbreakR.jaspar2022.besttfmodel.roc.prc.txt
#input for the script: -e: prefix of positive set and negative set,eg: evaldata_positive_data.txt,evaldata_negative_data.txt .The 'snp' and 'TF_SYMBOL' columns must be provided.
#-f: a file of merged experimental and predictive difference value(2 alleles of snp) of TF binding, the 'snp','TF_SYMBOL','model_name',predictive difference value of TF binding colums must be provided.
#-d: predictive difference value(2 alleles of snp) of TF binding

#6.2 HOCOMOCO v11
Rscript --vanilla auroc_auprc.R -e ../../snpdata/testdata/evaldata -f $datadir/motifbreakR.hocomocov11.merged.expe.pred.results.txt -m $datadir/evaldata_intermotifbreakR.hocomoco400tf.csv -d alleleDiff -o $datadir/motifbreakR.hocomocov11.tf.roc.prc.txt

# 7 regardless of pwm source,select TF's max AUROC and best models' merged experimental and predictive difference value(2 alleles off snp) of TF binding per TF,the data is used for computing correlation
Rscript --vanilla bestmotif.R --f1 $datadir/motifbreakR.jaspar2022.besttfmodel.roc.prc.txt --f2 $datadir/motifbreakR.hocomocov11.tf.roc.prc.txt --f3 $datadir/motifbreakR.jaspar2022.merged.expe.pred.results.txt --f4 $datadir/motifbreakR.hocomocov11.merged.expe.pred.results.txt --o1 $datadir/motifbreakR.alltf.bestmodel.roc.prc.txt --o2 $datadir/motifbreakR.alltf.merged.expe.pred.results.txt

