#TFs that motifbrekR can predict have intersection with selected TFs have at least 20 pbsnps.
#we only use common TFs for predicting snps' effect

#JASPAR 2022
modeltf<-'./JASPAR2022_691motif.csv'
evaldatatf<-'../snpdata/testdata/20pbsnptf.csv'
outfile<-'./testdata/evaldata_intermotifbreakR.jaspar587tf.csv'

df<-read.csv(modeltf,stringsAsFactors = F)
#intersection
tf<-read.csv(evaldatatf,stringsAsFactors=F)
df<-df[df$TF_SYMBOL %in% tf$TF_SYMBOL,]
write.csv(df,outfile,row.names=F,quote=F)

#HOCOMOCO v11
modeltf<-'./HOCOMOCOv11_400tf.csv'
evaldatatf<-'../snpdata/testdata/20pbsnptf.csv'
outfile<-'./testdata/evaldata_intermotifbreakR.hocomoco400tf.csv'

df<-read.csv(modeltf,stringsAsFactors = F)
#intersection
tf<-read.csv(evaldatatf,stringsAsFactors=F)
df<-df[df$TF_SYMBOL %in% tf$TF_SYMBOL,]
write.csv(df,outfile,row.names=F,quote=F)

