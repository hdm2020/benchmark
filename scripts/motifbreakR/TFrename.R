#motifbreakR: modify TF name as HGNC symbol

#JASPAR 2022
library(motifbreakR);library(MotifDb)
jaspar2022<-as.list(subset(MotifDb,dataSource=='jaspar2022' & organism=='Hsapiens'))
motifname<-data.frame(model_name=names(jaspar2022))#691
motifname$tf<-gsub('Hsapiens-jaspar2022-','',motifname$motif_name)
tfname<-function(xx){return(strsplit(xx,'-MA')[[1]][1])}
motifname$tf<-apply(data.frame(motifname$tf),1,tfname)
#modify TF name as standard TF name
library(org.Hs.eg.db)
df<-motifname
tfid<-mapIds(org.Hs.eg.db,df$tf,'ENTREZID','SYMBOL')
unique(df$tf[is.na(tfid)])#0
colnames(df)<-c('model_name','TF_SYMBOL')
write.csv(df,'JASPAR2022_691motif.csv',row.names=F,quote=F)#691tf,646tf,59tf-2tf intersection

#HOCOMOCO v11
hocomocov11<-as.list(subset(MotifDb,dataSource=='HOCOMOCOv11-core-A'|dataSource=='HOCOMOCOv11-core-B'|dataSource=='HOCOMOCOv11-core-C' & organism=='Hsapiens'))#400
motifname<-data.frame(model_name=names(hocomocov11))#400
tfinfo<-read.table('./HOCOMOCOv11_core_annotation_HUMAN_mono.tsv',header=T,sep='\t',stringsAsFactors=F)
tfinfo<-tfinfo[,c(1,2,7)]
tfname<-function(xx){return(strsplit(xx,'-')[[1]][5])}
motifname$Model<-apply(data.frame(motifname$motif_name),1,tfname)
motifname<-merge(motifname,tfinfo)
#modify TF name as standard TF name
motifname[motifname$Transcription.factor=='T','Transcription.factor']<-'TBXT'
colnames(motifname)[c(1,3)]<-c('model_name','TF_SYMBOL')
write.csv(motifname,'HOCOMOCOv11_400tf.csv',row.names=F,quote=F)

