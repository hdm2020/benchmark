#script for predicting snp's effect on TF binding using tRap

#construct a specific directory
datadir=testdata

#1 download tRap
wget http://trap.molgen.mpg.de/download/TRAP_R_package/tRap_0.5.tar.gz
#install tRap
R CMD INSTALL tRap_0.5.tar.gz


#2 learn the parameters of the generalized extreme value (GEV) distributions from an appropriate training set for using custom matrices
mkdir promoter_fa
cd promoter_fa
bash train_params.sh #note: the path new R package you should make a change


#3 modify names of TFs tRap can predict as their standard HGNC SYMBOL

#JASPAR 2022
wget -c https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt --no-check-certificate
grep '>' JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt > jaspar2022.name.txt
sed -i 's/>//g' jaspar2022.name.txt
sed -i '1i\motif_id\ttf' jaspar2022.name.txt

#HOCOMOCO v11
wget -c https://hocomoco11.autosome.org/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_HUMAN_mono_jaspar_format.txt --no-check-certificate
wget -c https://hocomoco11.autosome.org/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv --no-check-certificate

Rscript TFrename.R


#4 select some TFs to predict (TFs with >= 20 pbsnps for predicting)
mkdir $datadir
Rscript testdata.predictedTF.R #note: you need to change the code for your data


#5 model's input sequence: 40-bp
python generate_allelic_seqs.py -f ../../genome/hg19.fa -s ../../snpdata/testdata/testsnppos.tsv -o $datadir/testdata 2>$datadir/testdata.log
#input for the script: -f:reference genome; -s:a file of snps,a snp per line,eg:chr10_114258723_G_A


#6 model's prediction

#JASPAR 2022
mkdir $datadir/pwm_jaspar2022
mkdir $datadir/pwm_jaspar2022/results
nohup Rscript --vanilla tRap.predict.R -f $datadir/testdata -m $datadir/evaldata_intertRap.jaspar587tf.csv --pwmdb jaspar -o $datadir/pwm_jaspar2022/results/tRap.jaspar2022.txt &
#input for the script: -f:output of the script 'generate_allelic_seqs.py',prefix of the files including the sequences with ref or alt allele;
#-m: TF model you need to run,a TF model per line. The 'motif_name' column must be provided,eg:MA0069.1
#output for the script:a file including 11 columns: seq1, seq2, matrix, a1, a2, p1, p2, ratio, log.ratio, min.p, min.prod

#HOCOMOCO v11
mkdir $datadir/pwm_hocomocov11
mkdir $datadir/pwm_hocomocov11/results
nohup Rscript --vanilla tRap.predict.R -f $datadir/testdata -m $datadir/evaldata_intertRap.hocomoco401tf.csv --pwmdb hocomoco -o $datadir/pwm_hocomocov11/results/tRap.hocomocov11.txt &
#-m: TF model you need to run,a TF model per line. The 'motif_name' column must be provided,eg:ALX1_HUMAN.H11MO.0.B

#7 merge experimental and predictive difference value(2 alleles off snp) of TF binding
Rscript testdata.merge.R #note: you need to change the code for your data

#8 calculate AUROC,AUPRC of TFs

#8.1 JASPAR 2022
Rscript --vanilla auroc_auprc.R -e ../../snpdata/testdata/evaldata -f $datadir/tRap.jaspar2022.merged.expe.pred.results.txt -m $datadir/evaldata_intertRap.jaspar587tf.csv -d log.ratio -o $datadir/tRap.jaspar2022.tf.roc.prc.txt
mv besttfmodel.roc.prc.txt $datadir/tRap.jaspar2022.besttfmodel.roc.prc.txt
#input for the script: -e: prefix of positive set and negative set,eg: evaldata_positive_data.txt,evaldata_negative_data.txt .The 'snp' and 'TF_SYMBOL' columns must be provided.
#-f: a file of merged experimental and predictive difference value(2 alleles of snp) of TF binding, the 'snp','TF_SYMBOL','model_name',predictive difference value of TF binding colums must be provided.
#-d: predictive difference value(2 alleles of snp) of TF binding

#8.2 HOCOMOCO v11
Rscript --vanilla auroc_auprc.R -e ../../snpdata/testdata/evaldata -f $datadir/tRap.hocomocov11.merged.expe.pred.results.txt -m $datadir/evaldata_intertRap.hocomoco401tf.csv -d log.ratio -o $datadir/tRap.hocomocov11.tf.roc.prc.txt

#9 regardless of pwm source,select TF's max AUROC and best models' merged experimental and predictive difference value(2 alleles off snp) of TF binding per TF,the data is used for computing correlation
Rscript --vanilla bestmotif.R --f1 $datadir/tRap.jaspar2022.besttfmodel.roc.prc.txt --f2 $datadir/tRap.hocomocov11.tf.roc.prc.txt --f3 $datadir/tRap.jaspar2022.merged.expe.pred.results.txt --f4 $datadir/tRap.hocomocov11.merged.expe.pred.results.txt --o1 $datadir/tRap.alltf.bestmodel.roc.prc.txt --o2 $datadir/tRap.alltf.merged.expe.pred.results.txt

