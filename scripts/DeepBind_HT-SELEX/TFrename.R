#DeepBind_HT-SELEX:modify TF name as HGNC symbol
#awk -F '\t' '{print $1}' ~/.kipoi/model_names/DeepBind/model_names.tsv|grep 'Homo_sapiens/TF/'|grep 'SELEX' > hg.selex.378tf.model_names.txt
df<-read.table('hg.selex.378tf.models.txt',header=F,stringsAsFactors=F)
colnames(df)<-'model_name'
df$model_name<-paste0('DeepBind/',df$model_name)
tf<-function(xx){return(strsplit(xx,'_')[[1]][4])}
df$TF<-apply(df,1,tf)
#modify TF name as standard TF name
library(org.Hs.eg.db)
tfid<-mapIds(org.Hs.eg.db,df$TF,'ENTREZID','SYMBOL')
unique(df$TF[is.na(tfid)])
name1<-c("BHLHB2","BHLHB3","CART1","HINFP1","POU5F1P1","RAXL1","T","ZNF238","ZNF306","ZNF435")
name2<-c("BHLHE40","BHLHB3","ALX1","HINFP1","POU5F1B","RAX2","TBXT","ZBTB18","ZKSCAN3","ZSCAN16")#BHLHE41(BHLHB3) is existed.
df$TF_SYMBOL<-df$TF
for (i in 1:length(name1)){df[df$TF==name1[i],'TF_SYMBOL']<-name2[i]}
df<-df[df$model_name!='DeepBind/Homo_sapiens/TF/D00351.001_SELEX_EGR1',]#only EGR1 has 2 model_names,randomly select one
write.csv(df,'hg.selex.378tf.models.tf.csv',row.names=F,quote=F)

