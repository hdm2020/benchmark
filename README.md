# Comparative analysis of models in predicting the effects of SNPs on TF-DNA binding using large-scale in vitro and in vivo data

[TOC]

For users who want to run benchmarking pipeline or use these models, we provide step-by-step instructions for our comparative analysis with test data. Specifically, this instruction provides guidance and information on the following aspects:

1. Download test data and scripts
2. Install models, repository, packages and reference genome
3. Dependencies and requirements for environment and packages
4. Pre-process test data
5. Use provided scripts to predict SNP's effect on TF binding using different models and measure models' performance



## 1 Introduction

Alteration of transcription factor (TF) binding sites is one of the most important mechanisms for the functional impact of noncoding variations. Although many computational models are developed to predict the effects of noncoding variants on TF binding, there is a lack of systematic analysis to evaluate the predictive power of these methods. Given that, We evaluated 14 different models built on position weight matrices (PWMs), support vector machine (SVM), ordinary least squares (OLS) and deep neural networks (DNN) by large-scale in vitro (i.e. SNP-SELEX) and in vivo (i.e. allele-specific binding, ASB) TF binding data. 

4 PWM-based models (atSNP, motifbreakR, tRap and FABIAN-variant), 3 kmer/gkm-based machine learning methods (deltaSVM_HT-SELEX, deltaSVM_ ChIP-seq, QBiC-Pred), and 7 models based on DNN (DeepBind_HT-SELEX, DeepBind_ChIP-seq, DeepSEA, Beluga, DeepFun, Sei, Enformer) were incorporated in this study.



## 2 Download test data and scripts

You can directly download **all the files** in this github repository by the command below

```shell
git clone https://github.com/hdm2020/benchmark.git
```

Or you can use the following command

```shell
wget https://github.com/hdm2020/benchmark/archive/refs/heads/main.zip
unzip main.zip
```

If you are only interested in **certain directory**, you can download it using command:

```shell
# download snpdata/testdata
svn checkout https://github.com/hdm2020/benchmark/trunk/snpdata/testdata
# download scripts
svn checkout https://github.com/hdm2020/benchmark/trunk/scripts
```

If you are only interested in **certain file**, you can download it using command:

```shell
wget https://raw.githubusercontent.com/hdm2020/benchmark/main/<directoryName>/<FileName>
```



## 3 Software download and installation

Before you run the pipeline, please make sure that you have installed Python and R. We have already check Python 3.9.1 and R 4.2.0 is feasible for these models. We check the scripts in CentOS Linux release 7.9.2009 (Core). 

In order to run our scripts successfully, you need to download or install some models, repositories, packages and the reference genome. 

Table 1 below shows **the information for downloading each model**.

| Model             | Version            | Platform        | URL                                                          |
| ----------------- | ------------------ | --------------- | ------------------------------------------------------------ |
| atSNP             | atSNP_1.12.0       | R               | https://github.com/keleslab/atSNP                            |
| motifbreakR       | motifbreakR_2.10.0 | R               | https://github.com/Simon-Coetzee/MotifBreakR                 |
| tRap              | tRap_0.5           | R               | http://trap.molgen.mpg.de/                                   |
| deltaSVM_HT-SELEX |                    | C++;Perl;Python | https://github.com/ren-lab/deltaSVM                          |
| deltaSVM_ChIP-seq |                    | Python          | https://www.beerlab.org/deltasvm_models/                     |
| DeepBind_HT-SELEX |                    | Python          | http://kipoi.org/models/DeepBind/Homo_sapiens/TF/            |
| DeepBind_ChIP-seq |                    | Python          | http://kipoi.org/models/DeepBind/Homo_sapiens/TF/            |
| DeepSEA           |                    | Python          | http://kipoi.org/models/DeepSEA/predict/                     |
| Beluga            |                    | Python          | http://kipoi.org/models/DeepSEA/beluga/                      |
| Sei               |                    | Python          | https://github.com/FunctionLab/sei-framework/blob/main/      |
| Enformer          |                    | Python          | https://github.com/deepmind/deepmind-research/blob/master/enformer/ |

Table 2 below shows **the detailed requirements and feasible versions of the requirements for each model**.  Instructions 3.1-3.4 will guide you on how to download or install these and codes in the shell script **test.sh**  help you, too.

| Dependencies /Requirements | Type                        | Version                      | atSNP | motifbreakR | tRap | deltaSVM_HT-SELEX | deltaSVM_ChIP-seq | DeepBind_HT-SELEX | DeepBind_ChIP-seq | DeepSEA | Beluga | Sei  | Enformer |
| :------------------------: | --------------------------- | ---------------------------- | :---: | :---------: | :--: | :---------------: | :---------------: | :---------------: | :---------------: | :-----: | ------ | ---- | -------- |
|          hg19.fa           | reference genome            |                              |       |             |  ✔️   |         ✔️         |         ✔️         |         ✔️         |         ✔️         |    ✔️    | ✔️      |      | ✔️        |
|          optparse          | R package                   | optparse_1.7.3               |   ✔️   |      ✔️      |  ✔️   |                   |                   |                   |                   |         |        |      |          |
|           PRROC            | R package                   | PRROC_1.3.1                  |   ✔️   |      ✔️      |  ✔️   |         ✔️         |         ✔️         |         ✔️         |         ✔️         |    ✔️    | ✔️      | ✔️    | ✔️        |
|           kipoi            | application / python module | 0.8.0                        |       |             |      |                   |                   |         ✔️         |         ✔️         |    ✔️    | ✔️      |      |          |
|          kipoiseq          | python module               | 0.7.1 (0.5.2 for Enformer)   |       |             |      |                   |                   |         ✔️         |         ✔️         |    ✔️    | ✔️      |      | ✔️        |
|      libmamba-solver       | python module               |                              |       |             |      |                   |                   |         ✔️         |         ✔️         |    ✔️    | ✔️      |      |          |
|           numpy            | python module               | 1.20.3 (1.19.5 for Enformer) |       |             |      |                   |         ✔️         |         ✔️         |         ✔️         |    ✔️    | ✔️      | ✔️    | ✔️        |
|           pandas           | python module               | 1.3.1 (1.2.3 for Enformer)   |       |             |      |                   |         ✔️         |         ✔️         |         ✔️         |    ✔️    | ✔️      | ✔️    | ✔️        |
|           pysam            | python module               | 0.18.0                       |       |             |      |         ✔️         |         ✔️         |         ✔️         |         ✔️         |    ✔️    | ✔️      |      |          |
|          pytorch           | python module               | 1.9.0                        |       |             |      |                   |                   |                   |                   |         |        | ✔️    |          |
|           selene           | python module               | 0.5.0                        |       |             |      |                   |                   |                   |                   |         |        | ✔️    |          |
|           docopt           | python module               | 0.6.2                        |       |             |      |                   |                   |                   |                   |         |        | ✔️    |          |
|         dm-sonnet          | python module               | 2.0.1                        |       |             |      |                   |                   |                   |                   |         |        |      | ✔️        |
|       tensoflow-gpu        | python module               | 2.6.0                        |       |             |      |                   |                   |                   |                   |         |        |      | ✔️        |
|       tensorflow-hub       | python module               | 0.14.0                       |       |             |      |                   |                   |                   |                   |         |        |      | ✔️        |



### 	3.1 reference genome

When you **generate the models' input sequence**, you need a .fa file that contains the reference genome. In our scripts, we use **hg19.fa** which is downloaded from UCSC Genome Browser for the test data. Here is the command to download it:

```shell
mkdir genome
cd genome
wget -c ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
```

> In order to be consistent with the default path in the script, we create a new directory named **"genome"** which is at the same level as the directory "snpdata "and "scripts". Then we download hg19.fa in this directory. (You can also download hg38.fa through 'wget -c ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz' if needed)
>
> You can also download the reference genome file from other path or use the downloaded file. If you do so, remember to change the reference genome file path in the **4th step** of the shell script **test.sh**.



### 	3.2 R packages

Packages "**optparse**" and "**PRROC**" is required in R script **auroc_auprc.R**. If you haven't installed these packages , you can use the following command in R to install it:

```R
install.packages("optparse")
install.packages("PRROC")
```

Install  "**atSNP**"  package:

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("atSNP")

#dependencies
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19") #BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
```

Install  "**motifbreakR**"  package:

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("motifbreakR")

#dependencies
BiocManager::install("MotifDb")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19") #BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("BiocParallel")
```

Install  "**Biostrings**"  package for "**tRap**":

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
```



### 	3.3 kipoi, kipoiseq and libmamba-solver

If you want to run the scripts for **DeepBind**, **DeepSEA** and **Beluga** models, you need to install **kipoi** and **kipoiseq** first. 

Kipoi is a repository offering more than 2,000 individual trained models that cover key predictive tasks in genomics. Kipoiseq is used for standard set of data-loaders for training and making predictions for DNA sequence-based models. For more information about kipoi, you can visit kipoi and kipoiseq websites: https://github.com/kipoi/kipoi and https://github.com/kipoi/kipoiseq.

Here is an easy way to install Kipoi via pip

```
pip install kipoi
pip install kipoiseq
```

It is expected to run slower when installing models from kipoi without **libmamba solver**. So we recommend installing it using the following command: 

```shell
conda install conda-libmamba-solver
```



### 	3.4 Python modules

We can use pip to install modules "**numpy**", "**pandas**", "**pysam**" and so on.

```bash
pip install numpy
pip install pandas
pip install pysam
```



## 4 Pre-process test data

All steps to pre-process the data are in the **preprocess.R** file. You can directly run this R script using the command:

```shell
Rscript preprocess.R
```

> Before you run the R script  preprocess.R, you need to ensure the file GVATdb_novelbatch_100TF3000snp.csv and the file preprocess.R are in the **same directory**. Otherwise, you need to change the file path in the first command of preprocess.R :
>
> ```R
> df<-read.csv('/<DirectoryPath>/GVATdb_novelbatch_100TF3000snp.csv',stringsAsFactors=F)
> ```



## 5 Use provided scripts for the whole analysis

We provide the shell script named **"test.sh"** to make allelic TF binding predictions for any SNPs by 11 models (we use the web servers of "FABIAN-variant", "QBiC-Pred" and "DeepFun" for predicting SNPs' effect) and compute each TF's AUROC. 

The script **test.sh** of each model provide the pipeline of the whole analysis, mainly including: **1)** download and install the model. **2)** subset the list of available TFs a model can predict and modify TF names as their standard HGNC symbols. **3)** select some TFs you want to predict. **4)** SNPs as input: modify snps' format or extract DNA sequences with a specific length. **5)** predict SNPs' effect on TF binding. **6)** merge experimental and predictive allelic TF binding differences. **7)** compute the AUROC and AUPRC values of each TF and select the TF model with the highest AUROC if the TF has more than one models. **8)** best model's merged experimental and predictive allelic TF binding differences for each TF.  There are shell commands and annotations in the **test.sh** script. 

### 5.1 List of main scripts

| Author       | scripts                  | Language | purpose | input                                                        | output                                                       |
| ------------ | ------------------------ | -------- | ------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| Dongmei Han  | TFrename.R               | R        | 2       | see 5.2.1                                                    | see 5.2.1                                                    |
| Dongmei Han  | testdata.predictedTF.R   | R        | 3       | see 5.2.1                                                    | see 5.2.1                                                    |
| Yunjiang Qiu | generate_allelic_seqs.py | Python   | 4       | **-f**: a file with .fa suffix that contains reference genome sequence, eg: hg19.fa.<br />**-s**: a file of SNPs in which an SNP per line, eg: testsnppos.tsv<br />(Using deltaSVM_Chip-seq as an example) | **-o**: a directory path used to store output files, eg: testdata/results/<br />(Using deltaSVM_Chip-seq as an example) |
| Dongmei Han  | snp.R                    | R        | 4       | **-f**: a file of SNPs in which an SNP per line, eg: testsnppos.tsv<br /> | **-o**: a file of SNPs with new format.                      |
| Dongmei Han  | *.predict.R              | R        | 5       | **-f**: a file of SNPs in which an SNP per line, output of the script 'snp.R' for atSNP and motifbreakR; prefix of the files including the sequences with ref or alt allele, output of the script 'generate_allelic_seqs.py' for tRap.<br />**-m**: a file of TF models used for prediction. The 'motif_name' column must be provided.<br />**--pwmdb**: PWM database, "jaspar" or "hocomoco". <br />**-g**: genome name, eg: BSgenome.Hsapiens.UCSC.hg19.<br />(**-n**: number of cores that atSNP need, default=10.<br />**--method**: "default" or "ic" or "log" that motifbreakR provided, default="default".) | **-o**: a file of predictive SNP-TF binding.                 |
| Dongmei Han  | testdata.merge.R         | R        | 6       | see 5.2.1                                                    | see 5.2.1                                                    |
| Dongmei Han  | auroc_auprc.R            | R        | 7       | **-e**: prefix of positive and negative set, eg: evaldata (evaldata_positive_data.txt,evaldata_negative_data.txt). The 'snp' and 'TF_SYMBOL' columns must be provided.  <br />**-f**: a file of merged experimental and predictive allelic TF binding differences, eg:deltaSVM_ChIP-seq.merged.expe.pred.results.txt. The 'snp', 'TF_SYMBOL', 'model_name', and predictive difference value of TF binding columns must be provided.<br />**-m**: a file of TF models used for prediction, eg: evaldata_inter699e3tf.csv. The 'motif_name' column must be provided.<br />**-d**: measurement of predictive allelic TF binding difference, eg: deltaSVM.<br />(Using deltaSVM_Chip-seq as an example)<br /> | **-o**: a file of all TF models' AUROC and AUPRC values, eg: deltaSVM_ChIP-seq.tf.roc.prc.txt.<br />(Using deltaSVM_Chip-seq as an example, if a TF has multiple models, also output a file named 'besttfmodel.roc.prc.txt'  that includes the best models with the highest AUROC for each TF) |
| Dongmei Han  | bestmodel.R              | R        | 8       | **--f1**: a file of all TF models' best AUROC and AUPRC values, eg: deltaSVM_ChIP-seq.alltf.bestmodel.roc.prc.txt<br />**--f2**: a file of merged experimental and predictive allelic TF binding differences for each TF. eg: deltaSVM_ChIP-seq.merged.expe.pred.results.txt. | **--out**: a file of best models' merged experimental and predictive allelic TF binding differences for each TF, eg:deltaSVM_ChIP-seq.alltf.merged.expe.pred.results.txt. |
| Dongmei Han  | bestmotif.R              | R        | 8       | **--f1**: a file of all TF models' best AUROC and AUPRC values using JASPAR database, eg: motifbreakR.jaspar2022.besttfmodel.roc.prc.txt.<br />**--f2**: a file of all TF models' best AUROC and AUPRC values using HOCOMOCO v11 database, eg: motifbreakR.hocomocov11.tf.roc.prc.txt.<br />**--f3**: a file of merged experimental and predictive allelic TF binding differences for each TF using JASPAR database. eg: motifbreakR.jaspar2022.merged.expe.pred.results.txt.<br />**--f4**: a file of merged experimental and predictive allelic TF binding differences for each TF using HOCOMOCO v11 database. eg: motifbreakR.hocomocov11.merged.expe.pred.results.txt. | **--o1**: a file of all TF models' best AUROC and AUPRC values, regardless of pwm source. eg: motifbreakR.alltf.bestmodel.roc.prc.txt<br />**--o2**: a file of best models' merged experimental and predictive allelic TF binding differences for each TF, regardless of pwm source. eg: motifbreakR.alltf.merged.expe.pred.results.txt. |

Here we provide instructions in **5.2** about **running scripts directly** (you'd better run scripts step-by-step because some tasks run in the background using 'nohup...&' for saving time or you can modify the script). We also provide instructions in **5.3** about **running scripts step-by-step**, which is convenient to understand and resolve problems.

### 	5.2 Run directly

You can directly run the shell script after changing some codes in files using the command:

```shell
bash test.sh
```

If you want to use your own data,  you need to change the codes about data inputting for the  scripts below:

**testdata.predictedTF.R**, **testdata.merge.R** （for all models）; **train_params.sh** (for "tRap")



### 	5.3 An example of running step-by-step

#### 5.3.1 deltaSVM ChIP-seq

Here are the instructions about using  deltaSVM_ChIP-seq model for prediction according to the **test.sh** script.



The **1st step** is to download deltaSVM model and related script. The file size of `deltasvm_models_e3tf.tar.gz` is 9.2G, so it's maybe time-consuming.

```shell
wget -c http://www.beerlab.org/deltasvm_models/downloads/deltasvm_models_e3tf.tar.gz --no-check-certificate
mkdir e3tfmodels
tar -zxvf deltasvm_models_e3tf.tar.gz -C e3tfmodels/
rm deltasvm_models_e3tf.tar.gz
wget -c http://www.beerlab.org/deltasvm_models/downloads/score_snp_seq.py --no-check-certificate
```



The **2nd step** is to modify deltaSVM_ChIP-seq model's predictable TF names as standard HGNC symbols. Run the commands below directly to accomplish this step. 

input for **TFrename.R**: `tfe3.699.txt`, `699.models.txt`

output for **TFrename.R**: `tfe3.699model.csv`, a file that contains TF model information and includes 6 columns: 'model_name', 'source', 'TF', 'cell', 'TF_SYMBOL', 'weight'. 

'model_name': unique ID for each TF model. 'TF': original TF name. 'cell': cell type that each TF model trained on. 'TF_SYMBOL': standard HGNC symbol of TFs. 'weight': weight file of each TF.

```shell
wget -c http://www.beerlab.org/deltasvm_models/downloads/sample_list.txt --no-check-certificate
grep TF_E3 sample_list.txt > tfe3.699.txt
ls e3tfmodels/ > 699.models.txt
Rscript TFrename.R
```



The **3rd step** is to select some TFs to predict. We need to create a directory to store the output files generated using test data. Here we name this directory "testdata", you'd better create it at first.

input for **testdata.predictedTF.R**: `tfe3.699model.csv`, `20pbsnptf.csv`

`20pbsnptf.csv` is a file that contains TFs with >=20 positive SNPs. It includes 3 columns: 'TF', 'Count', 'TF_SYMBOL'. 'TF': original TF name. 'Count': count of positive SNPs per TF. 'TF_SYMBOL': standard HGNC symbol of TFs.

output for **testdata.predictedTF.R**: `evaldata_inter699e3tf.csv`, `predicted.tf.model.txt`. 

`evaldata_inter699e3tf.csv`  is a file that contains TF model information we select to predict. 'model_name': unique ID for each TF model. 'TF': original TF name. 'cell': cell type that each TF model trained on. 'TF_SYMBOL': standard HGNC symbol of TFs.

`predicted.tf.model.txt` contains the weight files of each TF model.

```shell
datadir=testdata 
mkdir $datadir
Rscript testdata.predictedTF.R
```

> 🚨 If you want to run models with your own evaluation data, you must replace  `20pbsnptf.csv` with your own file in **testdata.predictedTF.R**



The **4th step** is to generate the model's input sequence. Here we need a reference genome file, which can be downloaded from UCSC (see 3.1).

input: `hg19.fa`, `testsnppos.tsv` 

`testsnppos.tsv` is a file of SNPs, an SNP per line, eg:chr10_114258723_G_A

output: `testdata.ref.fa`, `testdata.ref.fa`, `testsnp.log`

`testdata.*.fa` contains genome sequences for each SNP.

```shell
python generate_allelic_seqs.py -f ../../genome/hg19.fa -s ../../snpdata/testdata/testsnppos.tsv -o $datadir/testdata 2>$datadir/testsnp.log
```



The **5th step** is to run model's prediction. We need to create a new directory to store the prediction results. 

input for **deltaSVM_ChIP-seq.predict.py**: `testdata.ref.fa`, `testdata.ref.fa`, `predicted.tf.model.txt`

output for **deltaSVM_ChIP-seq.predict.py**: predictive SNP-TF binding difference stored in files with names same as files in `predicted.tf.model.txt`.

```shell
mkdir $datadir/results
for weight in `cat testdata/predicted.tf.model.txt`
do
python deltaSVM_ChIP-seq.predict.py $datadir/testdata.ref.fa $datadir/testdata.alt.fa e3tfmodels/$weight $datadir/results/
done
```



The **6th step** is to merge experimental and predictive difference values (2 alleles of SNP) of TF binding. 

input: `evaldata_inter699e3tf.csv`, `GVAT_novelbatch_TFSYMBOL.txt`

`GVAT_novelbatch_TFSYMBOL.txt` is evaluation data that contains experimental SNP-TF binding differences.

output: `deltaSVM_ChIP-seq.merged.expe.pred.results.txt` , a file that contains merged experimental and predictive allelic TF binding differences.

```shell
Rscript testdata.merge.R
```

> 🚨 If you want to run models with your own evaluation data, you must replace `evaldata_inter699e3tf.csv` and  `GVAT_novelbatch_TFSYMBOL.txt`(evaluation data) with your own files in **testdata.merge.R**.



The **7th step** is to calculate AUROC, AUPRC values of TFs.

input: prefix of positive set and negative set, (`evaldata_positive_data.txt`,`evaldata_negative_data.txt`), `deltaSVM_ChIP-seq.merged.expe.pred.results.txt`, `evaldata_inter699e3tf.csv`

output: `deltaSVM_ChIP-seq.tf.roc.prc.txt`, `besttfmodel.roc.prc.txt` (output when TF have multiple models)

`deltaSVM_ChIP-seq.tf.roc.prc.txt` is a file that mainly includes each TF's AUROC, AUPRC, numbers of positive SNPs, model name, original name and standard HGNC symbol, cell types...

`besttfmodel.roc.prc.txt` is a file that contains each TF's best model with the highest AUROC value, information is the same as `deltaSVM_ChIP-seq.tf.roc.prc.txt` provides.

```shell
Rscript --vanilla auroc_auprc.R -e ../../snpdata/testdata/evaldata -f $datadir/deltaSVM_ChIP-seq.merged.expe.pred.results.txt -m $datadir/evaldata_inter699e3tf.csv -d deltaSVM -o $datadir/deltaSVM_ChIP-seq.tf.roc.prc.txt
mv besttfmodel.roc.prc.txt $datadir/deltaSVM_ChIP-seq.alltf.bestmodel.roc.prc.txt
```



The **8th step** is to output merged experimental and predictive difference value of SNP-TF binding per TF for the best model, which can be used for computing correlation. It's the last step in the **test.sh** script. 

input: `deltaSVM_ChIP-seq.alltf.bestmodel.roc.prc.txt`, `deltaSVM_ChIP-seq.merged.expe.pred.results.txt`

output: `deltaSVM_ChIP-seq.alltf.merged.expe.pred.results.txt`

```shell
Rscript --vanilla bestmodel.R --f1 $datadir/deltaSVM_ChIP-seq.alltf.bestmodel.roc.prc.txt --f2 $datadir/deltaSVM_ChIP-seq.merged.expe.pred.results.txt --out $datadir/deltaSVM_ChIP-seq.alltf.merged.expe.pred.results.txt
```



## 6 Reference

**Paper:**

atSNP: Zuo C, Shin S, Keleş S. atSNP: transcription factor binding affinity testing for regulatory SNP detection. *Bioinformatics*. 2015;31(20):3353-3355. https://doi.org/10.1093/bioinformatics/btv328

motifbreakR: Coetzee SG, Coetzee GA, Hazelett DJ. motifbreakR: an R/Bioconductor package for predicting variant effects at transcription factor binding sites. *Bioinformatics*. 2015;31(23):3847-3849. https://doi.org/10.1093/bioinformatics/btv470

tRap:Thomas-Chollier, M., Hufton, A., Heinig, M. *et al.* Transcription factor binding predictions using TRAP for the analysis of ChIP-seq data and regulatory SNPs. *Nat Protoc* **6**, 1860–1869 (2011). https://doi.org/10.1038/nprot.2011.409

FABIAN-variant: Steinhaus, R., Robinson, P. N., & Seelow, D. (2022). FABIAN-variant: predicting the effects of DNA variants on transcription factor binding. *Nucleic acids research*, *50*(W1), W322–W329. https://doi.org/10.1093/nar/gkac393

deltaSVM_HT-SELEX: Lee, D., Gorkin, D., Baker, M. *et al.* A method to predict the impact of regulatory variants from DNA sequence. *Nat Genet* **47**, 955–961 (2015). https://doi.org/10.1038/ng.3331

deltaSVM_ChIP-seq: Lee, D., Gorkin, D., Baker, M. *et al.* A method to predict the impact of regulatory variants from DNA sequence. *Nat Genet* **47**, 955–961 (2015). https://doi.org/10.1038/ng.3331

Shigaki D, Adato O, Adhikari AN, et al. Integration of multiple epigenomic marks improves prediction of variant impact in saturation mutagenesis reporter assay. *Hum Mutat*. 2019;40(9):1280-1291. https://doi.org/10.1002/humu.23797

QBiC-Pred: Martin, V., Zhao, J., Afek, A., Mielko, Z., & Gordân, R. (2019). QBiC-Pred: quantitative predictions of transcription factor binding changes due to sequence variants. *Nucleic acids research*, *47*(W1), W127–W135. https://doi.org/10.1093/nar/gkz363

DeepBind_HT-SELEX, DeepBind_ChIP-seq: Alipanahi, B., Delong, A., Weirauch, M. *et al.* Predicting the sequence specificities of DNA- and RNA-binding proteins by deep learning. *Nat Biotechnol* **33**, 831–838 (2015). https://doi.org/10.1038/nbt.3300

DeepSEA: Zhou, J., Troyanskaya, O. Predicting effects of noncoding variants with deep learning–based sequence model. *Nat Methods* **12**, 931–934 (2015). https://doi.org/10.1038/nmeth.3547

Beluga: Zhou, J., Theesfeld, C.L., Yao, K. *et al.* Deep learning sequence-based ab initio prediction of variant effects on expression and disease risk. *Nat Genet* **50**, 1171–1179 (2018). https://doi.org/10.1038/s41588-018-0160-6

DeepFun: Pei, G., Hu, R., Dai, Y., Manuel, A. M., Zhao, Z., & Jia, P. (2021). Predicting regulatory variants using a dense epigenomic mapped CNN model elucidated the molecular basis of trait-tissue associations. *Nucleic acids research*, *49*(1), 53–66. https://doi.org/10.1093/nar/gkaa1137

Sei: Chen, K.M., Wong, A.K., Troyanskaya, O.G. *et al.* A sequence-based global map of regulatory activity for deciphering human genetics. *Nat Genet* **54**, 940–949 (2022). https://doi.org/10.1038/s41588-022-01102-2

Enformer: Avsec, Ž., Agarwal, V., Visentin, D. *et al.* Effective gene expression prediction from sequence by integrating long-range interactions. *Nat Methods* **18**, 1196–1203 (2021). https://doi.org/10.1038/s41592-021-01252-x

**Tutorial/Code:**

atSNP: https://www.bioconductor.org/packages/release/bioc/vignettes/atSNP/inst/doc/atsnp-vignette.html

motifbreakR: https://bioconductor.org/packages/release/bioc/vignettes/motifbreakR/inst/doc/motifbreakR-vignette.html

tRap: http://trap.molgen.mpg.de/download/TRAP_R_package/tRap-tutorial.html

deltaSVM_HT-SELEX: https://github.com/ren-lab/deltaSVM/tree/master

deltaSVM_ChIP-seq: https://www.beerlab.org/deltasvm_models/

DeepBind_HT-SELEX: http://kipoi.org/models/DeepBind/Homo_sapiens/TF/

DeepBind_ChIP-seq: http://kipoi.org/models/DeepBind/Homo_sapiens/TF/

DeepSEA: http://kipoi.org/models/DeepSEA/predict/

Beluga: https://kipoi.org/models/DeepSEA/beluga/

Sei: https://github.com/FunctionLab/sei-framework/tree/main/

Enformer: https://github.com/deepmind/deepmind-research/blob/master/enformer/
