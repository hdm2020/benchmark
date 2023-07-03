#deltaSVM_HT-SELEX: modify TF name as HGNC symbol

#94 high confidence TFs
#ls gkmsvm_models/ > 94tf.model.txt
df<-read.table('94tf.model.txt',header=F,stringsAsFactors=F)
colnames(df)<-'weight'
tf<-function(xx){return(strsplit(xx,'_')[[1]][1])}
df$TF<-apply(df,1,tf)
#identify TF that is not standard name
library(org.Hs.eg.db)
tfid<-mapIds(org.Hs.eg.db,df$TF,'ENTREZID','SYMBOL')
unique(df$TF[is.na(tfid)])#'T'
df$TF_SYMBOL<-df$TF
df$TF_SYMBOL[df$TF_SYMBOL=='T']<-'TBXT'
df$model_name<-gsub('.model.txt','',df$weight,fixed=T)
write.csv(df,'94tf.model.csv',row.names=F,quote=F)

#all 533 TFs
#deltaSVM 533 tf model,remove 94 best model that have been selected.
#ls 533gkmsvm_models/ > 533gkmsvm_models.txt
df<-read.table('533gkmsvm_models.txt',header=F,stringsAsFactors=F)#5880 model
colnames(df)<-'weight'
tf<-function(xx){return(strsplit(xx,'_')[[1]][1])}
df$TF<-apply(df,1,tf)#533TF
library(org.Hs.eg.db)
tfid<-mapIds(org.Hs.eg.db,df$TF,'ENTREZID','SYMBOL')
unique(df$TF[is.na(tfid)])#"T"      "ZNF323" "ZSCAN5"
df$TF_SYMBOL<-df$TF
df$TF_SYMBOL[df$TF_SYMBOL=='T']<-'TBXT'
df$TF_SYMBOL[df$TF_SYMBOL=='ZNF323']<-'ZSCAN31'
df$TF_SYMBOL[df$TF_SYMBOL=='ZSCAN5']<-'ZSCAN5A'
df$model_name<-gsub('.model.txt','',df$weight,fixed=T)
write.csv(df,'533gkmsvm_models.tf.csv',row.names=F,quote=F)
#remove 94 selected best model
tf94<-read.csv('94tf.model.csv',stringsAsFactors=F)
df<-df[!(df$TF %in% tf94$TF),]#439tf,4556 model
write.csv(df,'439gkmsvm_models.tf.csv',row.names=F,quote=F)

