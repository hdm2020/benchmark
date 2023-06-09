#DeepSEA:modify TF name as HGNC symbol
#cp ~/.kipoi/models/DeepSEA/predictor.names ./
#grep 'DNase\|ac|\|me1|\|me2|\|me3|\|H2AZ' predictor.names > nonTF.list
df<-read.table('predictor.names',header=F,sep='|',stringsAsFactors=F)
colnames(df)<-c('cell_line','TF','treatment')
nontf<-read.table('nonTF.list',header = F,sep='|',stringsAsFactors = F)
df<-df[!(df$TF %in% nontf$V2),]
df$tf_coord<-as.numeric(row.names(df))
df$tfcoord_plus1<-df$tf_coord+1
df<-na.omit(df)#one line is NA
#modify TF name as standard TF name
library(org.Hs.eg.db)
tfid<-mapIds(org.Hs.eg.db,df$TF,'ENTREZID','SYMBOL')
unique(df$TF[is.na(tfid)])
df$TF_rename<-gsub('c-','',df$TF);df$TF_rename<-gsub('eGFP-','',df$TF_rename);df$TF_rename<-gsub('-','',df$TF_rename)
df$TF_rename<-toupper(df$TF_rename)
tfid<-mapIds(org.Hs.eg.db,df$TF_rename,'ENTREZID','SYMBOL')
name1<-c("JARID1A","POL2(B)","P300","PLU1","GABP","GR","NRSF","POL2","SIN3AK20","ERALPHA","PAX5C20","PAX5N19","POL24H8","PU.1","POL2(PHOSPHOS2)","NFKB","COREST","POL3","TBLR1","TR4","WHIP","KAP1","AP2ALPHA","AP2GAMMA","BAF155","BAF170","BRG1","HAE2F1","INI1","RPC155","SPT20","TFIIIC110","ERRA","GRP20","PGC1A","UBF")
name2<-c("KDM5A","POL2(B)","EP300","KDM5B","GABP","NR3C1","REST","POL2","SIN3AK20","ESR1","PAX5","PAX5","POL24H8","SPI1","POL2(PHOSPHOS2)","NFKB","RCOR1","POL3","TBL1XR1","NR2C2","WRNIP1","TRIM28","TFAP2A","TFAP2C","SMARCC1","SMARCC2","SMARCA4","HAE2F1","SMARCB1","POLR3A","SUPT20H","GTF3C2","ESRRA","GRP20","PPARGC1A","UBTF")
df$TF_SYMBOL<-df$TF_rename
for (i in 1:length(name1)){df[df$TF_rename==name1[i],'TF_SYMBOL']<-name2[i]}
df$model<-paste(df$cell_line,df$TF_SYMBOL,df$tf_coord,sep='_')
write.csv(df,'deepsea.689tf(1na).csv',row.names=F,quote=F)#163TF,689model(690 model,one is NA)

