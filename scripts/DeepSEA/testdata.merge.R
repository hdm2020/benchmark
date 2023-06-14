#merge experimental and predictive difference value(2 alleles off snp) of TF binding

modeltf<-'./testdata/evaldata_interdeepsea689tf.csv'
evaldatafile<-'../snpdata/testdata/GVAT_novelbatch_TFSYMBOL.txt'
predfile<-'./testdata/results/deepsea.txt'
outfile<-'./testdata/DeepSEA.merged.expe.pred.results.txt'

#read prediction result
model_tf<-read.csv(modeltf,stringsAsFactors=F)
result<-read.table(predfile,header=F,stringsAsFactors=F)
allresult<-data.frame()
for (i in 1:nrow(model_tf)){
  df<-result[,c(1,model_tf[i,'tfcoord_plus1'])]
  colnames(df)<-c('snp','delta_alt_ref')
  df$TF_SYMBOL<-model_tf[i,'TF_SYMBOL']
  df$model_name<-model_tf[i,'model_name']
  allresult<-rbind(allresult,df)
}


#read evaluation data
evaldata<-read.table(evaldatafile,sep='\t',header=T,stringsAsFactors = F)
evaldata<-evaldata[,c('snp','pbs','pval','TF_SYMBOL')]

#merge and output
df<-merge(evaldata,allresult,by=c('snp','TF_SYMBOL'))
write.table(df,outfile,col.names=T,row.names=F,quote=F,sep='\t')

