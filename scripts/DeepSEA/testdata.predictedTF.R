#TFs that DeepSEA can predict have intersection with selected TFs have at least 20 pbsnps.
#we only use common TFs for predicting snps' effect

modeltf<-'./deepsea.689tf(1na).csv'
evaldatatf<-'../../snpdata/testdata/20pbsnptf.csv'
outfile<-'./testdata/evaldata_interdeepsea689tf.csv'

df<-read.csv(modeltf,stringsAsFactors = F)
#intersection
tf<-read.csv(evaldatatf,stringsAsFactors=F)
df<-df[df$TF_SYMBOL %in% tf$TF_SYMBOL,]
write.csv(df,outfile,row.names=F,quote=F)

