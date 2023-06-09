#TFs that deltaSVM_HT-SELEX can predict have intersection with selected TFs have at least 20 pbsnps.
#we only use common TFs for predicting snps' effect

#94 high confidence TFs
modeltf<-'./94tf.model.csv'
evaldatatf<-'../snpdata/testdata/20pbsnptf.csv'
outfile<-'./testdata/evaldata_inter94tf.csv'

df<-read.csv(modeltf,stringsAsFactors = F)
#intersection
tf<-read.csv(evaldatatf,stringsAsFactors=F)
df<-df[df$TF_SYMBOL %in% tf$TF_SYMBOL,]
write.csv(df,outfile,row.names=F,quote=F)

df1<-df

#439 TFs
modeltf<-'./439gkmsvm_models.tf.csv'
evaldatatf<-'../snpdata/testdata/20pbsnptf.csv'
outfile1<-'./testdata/evaldata_inter439tf.csv'
outfile2<-'./testdata/predicted.tf.model.txt'
outfile3<-'./testdata/evaldata_inter533tf.csv'

df<-read.csv(modeltf,stringsAsFactors = F)
#intersection
tf<-read.csv(evaldatatf,stringsAsFactors=F)
df<-df[df$TF_SYMBOL %in% tf$TF_SYMBOL,]
write.csv(df,outfile1,row.names=F,quote=F)
write.table(df$model_name,outfile2,row.names=F,col.names=F,quote=F)

df2<-df

#all 533 TFs
alldf<-rbind(df1,df2)
write.csv(df,outfile3,row.names=F,quote=F)

