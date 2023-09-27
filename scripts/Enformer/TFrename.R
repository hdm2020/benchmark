#Enformer:modify TF name as HGNC symbol
library(openxlsx)
df<-read.xlsx('Enformer.targets.human.xlsx')
df<-df[df$assay_subtype=='ChIP-TF',c('index','identifier','description')]#2131
df$description<-gsub(',',';',df$description)
df$TF<-apply(data.frame(df$description),1,function(xx){return(strsplit(xx,':')[[1]][2])})#767

#modify TF name as standard TF name
#library(org.Hs.eg.db)
#tfid<-mapIds(org.Hs.eg.db,df$TF,'ENTREZID','SYMBOL')
#unique(df$TF[is.na(tfid)])
df$TF_rename<-gsub('eGFP-','',df$TF);df$TF_rename<-gsub('3xFLAG-','',df$TF_rename)
df<-df[df$TF!='.',]#these are not human,120
df<-df[df$TF!='RNAPII',]
df<-df[!(df$TF %in% c('H4K20me1','H2BK20ac','H4K8ac','H4K5ac','H2AK9ac','H2BK15ac','H2BK5ac','H4K91ac','H2BK12ac','H2BK120ac','H2AK5ac','H4K12ac')),]#
df[df$TF=='abcam','TF_rename']<-'BMAL1' #':' caused
df[df$TF=='active','TF_rename']<-'HIF1A' #':' caused #1914,TF_rename:685
#
name1<-c("H2AFZ","WHSC1","PTRF","C11orf30","ZNF645","CEBPb","GR","hBMAL1","hHIF1A")
name2<-c("H2AZ1","NSD2","CAVIN1","EMSY","CBLL2","CEBPB","NR3C1","BMAL1","HIF1A")
df$TF_SYMBOL<-df$TF_rename
for (i in 1:length(name1)){df[df$TF_rename==name1[i],'TF_SYMBOL']<-name2[i]}#1914,TF_SYMBOL:681TF(standard:678TF)

df$cell_line<-apply(data.frame(df$description),1,function(xx){return(strsplit(xx,':')[[1]][3])})


df$tf_coord<-df$index+1
df$tfcoord_plus1<-df$tf_coord+1

colnames(df)[2]<-'model_name'#unique id
df<-df[,2:9]
write.csv(df,'enformer.2131tf.csv',row.names=F,quote=F)#678TF,1914model

