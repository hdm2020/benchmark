#select every tf's best motif from JASPAR2022 and HOCOMOCO V11 database

#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option("--f1", type="character", default=NULL,
              help="auroc file name(jaspar)", metavar="character"),
  make_option("--f2", type="character", default=NULL,
              help="auroc file name(hocomoco)", metavar="character"),
  make_option("--f3", type="character", default=NULL,
              help="the file including merged experimental and predictive difference value of TF binding(jaspar)", metavar="character"),
  make_option("--f4", type="character", default=NULL,
              help="the file including merged experimental and predictive difference value of TF binding(hocomoco)", metavar="character"),
  make_option("--o1", type="character", default=NULL,
              help="output file name:auroc file name(jaspar & hocomoco)", metavar="character"),
  make_option("--o2", type="character", default=NULL,
              help="output file name:merged experimental and predictive difference value of TF binding(jaspar & hocomoco)", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#1 AUROC
#JASPAR 2022
jaspar<-read.table(opt$f1,header=T,sep='\t',stringsAsFactors=F)
#HOCOMOCO V11
hocomoco<-read.table(opt$f2,header=T,sep='\t',stringsAsFactors=F)

sharetf<-intersect(jaspar$TF_SYMBOL,hocomoco$TF_SYMBOL)
df<-rbind(jaspar[!(jaspar$TF_SYMBOL %in% sharetf),],hocomoco[!(hocomoco$TF_SYMBOL %in% sharetf),])
for (tf in sharetf){
  if (jaspar[jaspar$TF_SYMBOL==tf,'auroc']>=hocomoco[hocomoco$TF==tf,'auroc']){
    df<-rbind(df,jaspar[jaspar$TF_SYMBOL==tf,])}
  else{
    df<-rbind(df,hocomoco[hocomoco$TF_SYMBOL==tf,])}}
write.table(df,opt$o1,row.names=F,col.names=T,sep='\t',quote=F)

#2 merged experimental and predictive difference value of TF binding
#JASPAR 2022
jaspar<-read.table(opt$f3,header=T,sep='\t',stringsAsFactors=F)
#HOCOMOCO V11
hocomoco<-read.table(opt$f4,header=T,sep='\t',stringsAsFactors=F)

result<-rbind(jaspar[jaspar$model_name %in% df$model_name,],hocomoco[hocomoco$model_name %in% df$model_name,])
write.table(result,opt$o2,row.names=F,col.names=T,sep='\t',quote=F)

