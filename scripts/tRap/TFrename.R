#tRap: modify TF name as HGNC symbol

#JASPAR 2022
#download pfm matrix of human tfs from jaspar 2022.
#wget -c https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt --no-check-certificate
#grep '>' JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt > jaspar2022.name.txt
#sed -i 's/>//g' jaspar2022.name.txt
#sed -i '1i\motif_id\ttf' jaspar2022.name.txt
###subset human motif according to the strings and search in JASPAR database
motif<-read.table('jaspar2022.name.txt',header=T,sep='\t',stringsAsFactors=F)#841 motif
motif$TF<-toupper(motif$tf)
df<-motif[motif$tf==motif$TF,]#695 motif,650tf,60 are 2tf-intersection
df<-df[!(df$motif_id %in% c('MA0108.2','MA0259.1','MA0619.1','MA1540.2')),]#MA0108.2 TBP(no species)；MA0259.1 ARNT::HIF1A(multi-species)；MA0619.1 LIN54(Gallus gallus)；MA1540.2 NR5A1(Mus musculus)
#modify TF name as standard HGNC symbol
library(org.Hs.eg.db)
tfid<-mapIds(org.Hs.eg.db,df$tf,'ENTREZID','SYMBOL')
unique(df$tf[is.na(tfid)])#0
colnames(df)<-c('model_name','TF_SYMBOL')
write.csv(df[,1:2],'JASPAR2022_human587tf.csv',row.names=F,quote=F)

#HOCOMOCO v11
#download jaspar format
#wget -c https://hocomoco11.autosome.org/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_HUMAN_mono_jaspar_format.txt --no-check-certificate
#add a blank line at the end,because there is a warning during the process of loading the motif library,but it can't influence result
#download annotation information
#wget -c https://hocomoco11.autosome.org/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv --no-check-certificate
df<-read.table('HOCOMOCOv11_core_annotation_HUMAN_mono.tsv',header=T,sep='\t',stringsAsFactors=F)
df<-df[,c(1,2)]
#modify TF name as standard TF name
library(org.Hs.eg.db)
tfid<-mapIds(org.Hs.eg.db,df$Transcription.factor,'ENTREZID','SYMBOL')
unique(df$Transcription.factor[is.na(tfid)])#'T'
df[df$Transcription.factor=='T','Transcription.factor']<-'TBXT'
colnames(df)[1:2]<-c('model_name','TF_SYMBOL')
write.csv(df,'HOCOMOCOv11_401tf.csv',row.names=F,quote=F)

