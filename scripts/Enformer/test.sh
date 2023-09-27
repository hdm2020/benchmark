#script for predicting snp's effect on TF binding using Enformer

#construct a specific directory
datadir=testdata

#1 download Enformer
wget -c https://raw.githubusercontent.com/deepmind/deepmind-research/master/enformer/attention_module.py
wget -c https://raw.githubusercontent.com/deepmind/deepmind-research/master/enformer/enformer.py
wget -c https://raw.githubusercontent.com/deepmind/deepmind-research/master/enformer/enformer_test.py

#GPU version
python3 -m venv enformer_gpu_venv
source enformer_gpu_venv/bin/activate

#install dependcies
pip install dm-sonnet #(2.0.1)
pip install numpy==1.19.5
pip install pandas==1.2.3
pip install kipoiseq==0.5.2
pip install tensorflow-gpu==2.6.0 -i https://pypi.tuna.tsinghua.edu.cn/simple
pip install tensorflow-hub #(0.14.0)
#pip uninstall tensorflow-estimator (2.13.0)
#pip install tensorflow-estimator==2.6.0
#pip uninstall protobuf(4.23.4)
#pip install protobuf==3.20.0
#pip uninstall keras(2.13.1)
#pip install keras==2.6.0
conda install cudnn #Could not load dynamic library 'libcudnn.so.8'

#if you want to use jupyter,install ipykernel package
#pip install ipykernel
#add virtual environment into jupyter and set name
#python -m ipykernel install --user --name enformer_gpu

#download module 
mkdir moduleE
curl -L "https://hub.tensorflow.google.cn/deepmind/enformer/1?tf-hub-format=compressed" -k | tar -zxvC ./moduleE #-k/--insecure:Tell libcurl to not verify the peer; we download it locally beacause other URL can't be downloaded successfully.


#2 modify names of TFs Enformer can predict as their standard HGNC SYMBOL
#detailed Human feature information are downloaded from original paper-Supplementary Tables 2,and we keep first 11 columns of the file as Enformer.targets.human.xlsx.
#wget -c https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-021-01252-x/MediaObjects/41592_2021_1252_MOESM3_ESM.xlsx
Rscript TFrename.R


#3 select some TFs to predict (TFs with >= 20 pbsnps for predicting)
mkdir $datadir
Rscript testdata.predictedTF.R #note: you need to change the code for your data


#4 model's input
#a file of snps,a snp per line,eg:chr10_114258723_G_A


#5 model's prediction
mkdir $datadir/results

nohup python enformer.predict.py -s ../snpdata/testdata/testsnppos.tsv -g ../genome/hg19 -o $datadir/results/Enformer.txt > $datadir/testdata.log 2>&1 &

#6 merge experimental and predictive difference value(2 alleles of snp) of TF binding
Rscript testdata.merge.R #note: you need to change the code for your data


#7 calculate AUROC,AUPRC of TFs

Rscript --vanilla auroc_auprc.R -e ../snpdata/testdata/evaldata -f $datadir/Enformer.merged.expe.pred.results.txt -m $datadir/evaldata_interenformer2131tf.csv -d delta_alt_ref -o $datadir/Enformer.tf.roc.prc.txt
mv besttfmodel.roc.prc.txt $datadir/Enformer.alltf.bestmodel.roc.prc.txt
#input for the script: -e: prefix of positive set and negative set,eg: evaldata_positive_data.txt,evaldata_negative_data.txt .The 'snp' and 'TF_SYMBOL' columns must be provided.
#-f: a file of merged experimental and predictive difference value(2 alleles of snp) of TF binding, the 'snp','TF_SYMBOL','model_name',predictive difference value of TF binding colums must be provided.
#-d: predictive difference value(2 alleles of snp) of TF binding


#8 best models' merged experimental and predictive difference value(2 alleles off snp) of TF binding per TF,the data is used for computing correlation

Rscript --vanilla bestmodel.R --f1 $datadir/Enformer.alltf.bestmodel.roc.prc.txt --f2 $datadir/Enformer.merged.expe.pred.results.txt --out $datadir/Enformer.alltf.merged.expe.pred.results.txt

