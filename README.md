# Systematic comparative analysis of models on predicting the effects of SNPs on TF-DNA binding using large-scale in vitro and in vivo data
# Introduction
Although many computational models are developed to predict the effects of noncoding variants on TF binding, there is a lack of systematic analysis to evaluate the predictive power of these methods. Given that, We evaluated 10 different models built on position weight matrices (PWMs), support vector machine (SVM), ordinary least squares (OLS) and deep neural networks (DNN) by large-scale in vitro (SNP-SELEX) and in vivo (allele-specific binding, ASB) TF binding data. There are 3 PWM-based models (atSNP, motifbreakR, and tRap), 3 kmer/gkm-based machine learning methods (deltaSVM_HT-SELEX, deltaSVM_ ChIP-seq, QBiC-Pred), and 4 models based on DNN (DeepBind_HT-SELEX, DeepBind_ChIP-seq, DeepSEA, DeepFun). 
# Dependency
pysam is required to run deltaSVM_HT-SELEX model, kipoi is required to run DeepBind and DeepSEA. You can install them using conda or pip.
# Quick Start
Here we provide scripts to make allelic TF binding predictions for any SNPs by 8 models (we use web server of QBiC-Pred and DeepFun for predicting SNPs' effect) and compute each TF's AUROC. The script test.sh of each model provide the pipeline, mainly including 1) subset the list of available TFs a model can predict and modify TF names as their standard HGNC symbols. 2) select some TFs you want to predict. 3) SNPs as input: modify snps' format or extract DNA sequences with a specific length. 4) predict SNPs' effect on TF binding. 5) merge experimental and predictive allelic TF binding difference. 6) compute AUROC and AUPRC values of each TF. 
Before the prediction, you should make a preprocess for your benchmarking data to obtain TFs with more than 20 positive samples, SNPs and positive set and negative set. In addition, you can download reference genome from the UCSC Genome Browser.
NOTE: some codes may need to change according to your data, we make a hint in the script test.sh.
