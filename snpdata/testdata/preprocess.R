#subset TFs have at least 20 pbsnp for sufficient statistical power;
#subset TFs and snps of snp-tf pairs
#positive samples and negative samples for calculating AUROC,AUPRC

#1 subset TFs have at least 20 pbsnp(potential binding snp);and modify TF name as standard HGNC symbol(GeneCards)

df<-read.csv('GVATdb_novelbatch_100TF3000snp.csv',stringsAsFactors=F)

#modify TF name (we check these TFs with >=1 pbsnp)
xx<-data.frame(table(df[df$pval<0.01,'TF']))
xx$Var1<-as.character(xx$Var1)

library(org.Hs.eg.db)
tfid<-mapIds(org.Hs.eg.db,xx$Var1,'ENTREZID','SYMBOL');unique(xx$Var1[is.na(tfid)])
tf<-c("BHLHB2","ENSG00000206218","FOXO3A")#3 TFs don't have stantard symbols. 
tf_symbol<-c("BHLHE40","ENSG00000206218","FOXO3")
tf_tfsymbol<-data.frame(cbind(tf,tf_symbol))
#write.table(tf_tfsymbol,'tf_symbol.txt',row.names = F,col.names = T,quote = F,sep = '\t')
colnames(xx)<-c('TF','Count');xx$TF_SYMBOL<-xx$TF
for (i in 1: nrow(tf_tfsymbol)){xx$TF_SYMBOL[xx$TF==tf_tfsymbol[i,'tf']]<-tf_tfsymbol[i,'tf_symbol']}

#TFs with >=20pbsnp 
tf20<-subset(xx,xx$Count>=20)
write.csv(tf20,'20pbsnptf.csv',row.names = F,quote = F)


#2 add TF_SYMBOL column in the data

df$TF_SYMBOL<-df$TF
for (i in 1: nrow(tf_tfsymbol)){df$TF_SYMBOL[df$TF==tf_tfsymbol[i,'tf']]<-tf_tfsymbol[i,'tf_symbol']}
write.table(df[,c(1,2,5:7)],'GVAT_novelbatch_TFSYMBOL.txt',row.names = F,col.names = T,quote = F,sep = '\t')


#3 positive samples and negative samples for calculating AUROC,AUPRC
poss<-df[df$pval<0.01,]
write.table(poss[,c(1,2,5:7)],'evaldata_positive_data.txt',row.names = F,col.names = T,quote = F,sep = '\t')
negs<-df[df$pval>0.5,]
write.table(negs[,c(1,2,5:7)],'evaldata_negative_data.txt',row.names = F,col.names = T,quote = F,sep = '\t')


#4 snps
write.table(unique(df$snp),'testsnppos.tsv',row.names = F,col.names = F,quote = F,sep = '\t')
#awk -F '\t' '{print $2}' GVAT_novelbatch_TFSYMBOL.txt|sort|uniq > testsnppos.tsv
#sed -i '3001d' testsnppos.tsv

