#deltaSVM_ChIP-seq: modify TF name as HGNC symbol
#wget -c http://www.beerlab.org/deltasvm_models/downloads/sample_list.txt --no-check-certificate
#grep TF_E3 sample_list.txt >tfe3.699.txt
#ls e3tfmodels/ > 699.models.txt
df<-read.delim('tfe3.699.txt',header=F,stringsAsFactors=F)
df<-df[,c(1,3)];colnames(df)<-c('model_name','source')
for (i in 1:699){
df$TF[i]<-strsplit(df$source[i],' ')[[1]][1]
df$cell[i]<-strsplit(df$source[i],' ')[[1]][5]
}
df$TF<-gsub('eGFP-','',df$TF)
df$TF<-gsub('3xFLAG-','',df$TF)
#modify TF name as standard TF name
#library(org.Hs.eg.db)
#tfid<-mapIds(org.Hs.eg.db,df$TF,'ENTREZID','SYMBOL')
#unique(df$TF[is.na(tfid)])#"C11orf30"
df$TF_SYMBOL<-df$TF
df[df$TF=='C11orf30','TF_SYMBOL']<-'EMSY'
#merge weight file name
models<-read.table('699.models.txt',header=F,stringsAsFactors=F)
colnames(models)<-'weight'
for (i in 1:699){
models$model_name[i]<-strsplit(models$weight[i],'_hg38')[[1]][1]
}
df<-merge(df,models)
write.csv(df,'tfe3.699model.csv',row.names=F,quote=F)#465tf,699 model

