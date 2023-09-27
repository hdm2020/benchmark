#Sei:modify TF name as HGNC symbol

df<-read.delim('model/target.names',header=F,sep='\n',stringsAsFactors=F)
df$V1<-gsub(',_','/',fixed=T,df$V1)
df$V1<-gsub(',','/',fixed=T,df$V1)
df$cell_line<-apply(data.frame(df$V1),1,function(xx){return(strsplit(xx,' | ',fixed=T)[[1]][1])})
df$TF<-apply(data.frame(df$V1),1,function(xx){return(strsplit(xx,' | ',fixed=T)[[1]][2])})
df$source<-apply(data.frame(df$V1),1,function(xx){return(strsplit(xx,' | ',fixed=T)[[1]][3])})

#position index of histone marks('histone_inds.npy',0-based)
#import numpy as np
#import pandas as pd
#histone_inds = np.load('model/histone_inds.npy')
#df=pd.DataFrame(histone_inds)
#df.to_csv('histone_inds.txt',index=False,header=False,sep='\t')
histone_inds<-read.table('histone_inds.txt',header=F,stringsAsFactors=F)
histone_inds<-histone_inds+1
df<-df[setdiff(1:nrow(df),histone_inds$V1),]

#chromatin accessibility
#grep 'DNase\|ATAC-seq' model/target.names > chrom_access.txt
acce<-read.delim('chrom_access.txt',header=F,sep='\n',stringsAsFactors=F)
acce$feature<-apply(data.frame(acce$V1),1,function(xx){return(strsplit(xx,' | ',fixed=T)[[1]][2])})
df<-df[!(df$TF %in% unique(acce$feature)),]

df$tf_coord<-as.numeric(row.names(df))
write.table(df[,2:4],'Sei.TF.binding.txt',row.names=F,col.names=T,sep='\t',quote=F)

#modify TF name as standard TF name
library(org.Hs.eg.db)
df$TF_rename<-gsub('c-','',df$TF);df$TF_rename<-gsub('eGFP-','',df$TF_rename);df$TF_rename<-gsub('-','',df$TF_rename)
df$TF_rename<-toupper(df$TF_rename)
tfid<-mapIds(org.Hs.eg.db,df$TF_rename,'ENTREZID','SYMBOL')
unique(df$TF_rename[is.na(tfid)])#74
name1<-c('C17ORF49','CCDC101','T','PR','NKX21','NKX31','NKX22','C17ORF96','STAT5','WHSC1','CPSF3L','MRE11A','PTRF','MBD1_ISOFORM2','MBD1_ISOFORM1','C11ORF30','ZNF788','GUCY1B3','MGEA5','TAZ','MKL2','MKL1','FAM208A','GR','NRSF','P300','ERALPHA','COREST','PAX5C20','PAX5N19','PU.1','TBLR1','TR4','WHIP','JARID1A','KAP1','AP2ALPHA','AP2GAMMA','BAF155','BAF170','BRG1','INI1','RPC155','SPT20','TFIIIC110','ERRA','PGC1A','SREBP1','PLU1','UBF')
name2<-c('C17orf49','SGF29','TBXT','PGR','NKX2-1','NKX3-1','NKX2-2','EPOP','STAT5A','NSD2','INTS11','MRE11','CAVIN1','MBD1','MBD1','EMSY','ZNF788P','GUCY1B1','OGA','TAFAZZIN','MRTFB','MRTFA','TASOR','NR3C1','REST','EP300','ESR1','RCOR1','PAX5','PAX5','SPI1','TBL1XR1','NR2C2','WRNIP1','KDM5A','TRIM28','TFAP2A','TFAP2C','SMARCC1','SMARCC2','SMARCA4','SMARCB1','POLR3A','SUPT20H','GTF3C2','ESRRA','PPARGC1A','SREBF1','KDM5B','UBTF')
df$TF_SYMBOL<-df$TF_rename
for (i in 1:length(name1)){df[df$TF_rename==name1[i],'TF_SYMBOL']<-name2[i]}
df$model_name<-paste(df$cell_line,df$TF_SYMBOL,df$tf_coord,sep='_')
write.csv(df[,2:8],'Sei.9315tfmodel.csv',row.names=F,quote=F)#1089TF,1030TF SYMBOL(1006 HGNC symbol),9315 model


