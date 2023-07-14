#DeepBind_ChIP-seq:modify TF name as HGNC symbol
#awk -F '\t' '{print $1}' ~/.kipoi/models/DeepBind/models.tsv|grep 'Homo_sapiens/TF/'|grep 'ChIP-seq' > hg.chipseq.137tf.models.txt
df<-read.table('hg.chipseq.137tf.models.txt',header=F,stringsAsFactors=F)
colnames(df)<-'model_name'
df$model_name<-paste0('DeepBind/',df$model_name)
tf<-function(xx){return(strsplit(xx,'_')[[1]][4])}
df$TF<-apply(df,1,tf)#136tf,137 models
#modify TF name as standard TF name
#library(org.Hs.eg.db)
#tfid<-mapIds(org.Hs.eg.db,df$TF,'ENTREZID','SYMBOL')
#unique(df$TF[is.na(tfid)])
name1<-c("FAM48A","KAP1","RPC155","SIN3AK20")
name2<-c("SUPT20H","KAP1","POLR3A","SIN3AK20")#TRIM28(KAP1) is existed.
df$TF_SYMBOL<-df$TF
for (i in 1:length(name1)){df[df$TF==name1[i],'TF_SYMBOL']<-name2[i]}
df<-df[df$model_name!='DeepBind/Homo_sapiens/TF/D00755.005_ChIP-seq_EP300',]#only EP300 has 2 models,randomly select one
write.csv(df,'hg.chipseq.137tf.models.tf.csv',row.names=F,quote=F)
