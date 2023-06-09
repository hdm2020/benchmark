#TFs that deltaSVM_ChIP-seq can predict have intersection with selected TFs have at least 20 pbsnps.
#we only use common TFs for predicting snps' effect

modeltf<-'./tfe3.699model.csv'
evaldatatf<-'../snpdata/testdata/20pbsnptf.csv'
outfile1<-'./testdata/evaldata_inter699e3tf.csv'
outfile2<-'./testdata/predicted.tf.model.txt'

df<-read.csv(modeltf,stringsAsFactors = F)
#intersection
tf<-read.csv(evaldatatf,stringsAsFactors=F)
df<-df[df$TF_SYMBOL %in% tf$TF_SYMBOL,]
write.csv(df,outfile1,row.names=F,quote=F)
write.table(df$weight,outfile2,row.names=F,col.names=F,quote=F)

