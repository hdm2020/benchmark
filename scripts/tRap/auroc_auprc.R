#the script for calculating every TF's AUROC,AUPRC

#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-e", "--evaldata"), type="character", default=NULL,
              help="evaluation data (the files including positive samples and negative samples)", metavar="character"),
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="the file including merged experimental and predictive difference value of TF binding", metavar="character"),
  make_option(c("-m", "--model"), type="character", default=NULL,
              help="the file including model file name:selected TF model for predicting", metavar="character"),
  make_option(c("-d", "--deltascore"), type="character", default=NULL,
              help="predictive difference value of TF binding", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$evaldata)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
if (is.null(opt$model)){
  print_help(opt_parser)
  stop("You need to provide TF model list you have predicted.n", call.=FALSE)
}
if (is.null(opt$deltascore)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (predictive difference value of TF binding).n", call.=FALSE)
}


library(PRROC)
#define function
rocprc<-function(inputfile1,inputfile2,modeltf,deltascore,outputfile){
  #positive samples
  poss<-read.table(paste0(inputfile1,'_positive_data.txt'),header = T,stringsAsFactors = F,sep='\t')[,c('snp','TF_SYMBOL')]
  poss$clables<-1
  
  #negative samples
  negs<-read.table(paste0(inputfile1,'_negative_data.txt'),header = T,stringsAsFactors = F,sep='\t')[,c('snp','TF_SYMBOL')]
  negs$clables<-0
  
  result<-rbind(poss,negs)
  
  #merged experimental and predictive difference value of TF binding
  exp_pre<-read.table(inputfile2,header = T,stringsAsFactors = F,sep='\t')
  result<-merge(result,exp_pre)
   
  #absolute values of predictive difference values
  result<-result[!is.na(result[,deltascore]),]#2 conditions:motif's parameters is missed or p value=0---note:some motifs will be removed
  result[,deltascore]<-abs(result[,deltascore])
  #print(head(result))  
  
  #model_tf file
  model_tf<-read.csv(modeltf,stringsAsFactors = F)
  models<-unique(model_tf$model_name)
  n<-length(models)
  
  #AUROC,AUPRC
  auroc_auprc<-data.frame(auroc=numeric(n),auprc=numeric(n),num=numeric(n),row.names = models)
  for (model in models){
    df<-subset(result,result$model_name==model)
    n_table<-dim(as.data.frame(table(df$clables)))[1] #check if positive or negative samples of the TF are existed
    if (n_table==2){
      set.seed(1)
      m<-table(df$clables)[['1']] #the count of positive samples
      neg<-data.frame(neg_delta=sample(df[df$clables==0,deltascore],m,replace=F)) #randomly sample negative samples (same number as positive samples) 
      neg$clables_neg<-0
      neg<-cbind(neg,df[df$clables==1,c(deltascore,'clables')])

      auroc<-roc.curve(scores.class0 = c(neg$neg_delta,neg[,deltascore]), weights.class0 = c(neg$clables_neg,neg$clables))$auc
      auprc<-pr.curve(scores.class0 = c(neg$neg_delta,neg[,deltascore]), weights.class0 = c(neg$clables_neg,neg$clables))$auc.integral
      auroc_auprc[model,]<-c(auroc,auprc,m)
    }}
  auroc_auprc$model_name<-row.names(auroc_auprc)
  auroc_auprc<-merge(auroc_auprc,model_tf)
  
  #save results
  write.table(auroc_auprc,outputfile,col.names=T,row.names=F,quote=F,sep = '\t')
  
  #if a TF has more than one models,we select the model that has maximal AUROC
  n1=length(unique(auroc_auprc$model_name))
  n2=length(unique(auroc_auprc$TF_SYMBOL))
  if (n1!=n2){
    tf_num<-data.frame(table(auroc_auprc$TF_SYMBOL))
    tfrepeat<-tf_num[tf_num$Freq>1,]
    xx<-auroc_auprc[auroc_auprc$TF_SYMBOL %in% tfrepeat$Var1,]
    auroc_auprc<-auroc_auprc[!(auroc_auprc$TF_SYMBOL %in% tfrepeat$Var1),]
    for (tf in tfrepeat$Var1){
      tfs<-xx[xx$TF_SYMBOL==tf,]
      tf_maxauc<-tfs[tfs$auroc==max(tfs$auroc),][1,]
      auroc_auprc<-rbind(auroc_auprc,tf_maxauc)
    }
  write.table(auroc_auprc,'besttfmodel.roc.prc.txt',row.names=F,col.names=T,quote=F,sep='\t')
}}

rocprc(opt$evaldata,opt$file,opt$model,opt$deltascore,opt$out)

