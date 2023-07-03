#merge experimental and predictive difference value(2 alleles of snp) of TF binding

modeltf<-'./testdata/evaldata_interselex377tf.csv'
evaldatafile<-'../../snpdata/testdata/GVAT_novelbatch_TFSYMBOL.txt'
outfile<-'./testdata/DeepBind_HT-SELEX.merged.expe.pred.results.txt'

#read prediction result
model_tf<-read.csv(modeltf,stringsAsFactors=F)
result<-data.frame()
for (i in 1:nrow(model_tf)){
  df<-read.table(paste0('testdata/results/',model_tf[i,'TF']),header=T,stringsAsFactors=F,sep='\t')
  df$TF_SYMBOL<-model_tf[i,'TF_SYMBOL']
  df$model_name<-model_tf[i,'model_name']
  result<-rbind(result,df)
}

#read evaluation data
evaldata<-read.table(evaldatafile,sep='\t',header=T,stringsAsFactors = F)
evaldata<-evaldata[,c('snp','pbs','pval','TF_SYMBOL')]

#merge and output
df<-merge(evaldata,result,by=c('snp','TF_SYMBOL'))
write.table(df,outfile,col.names=T,row.names=F,quote=F,sep='\t')

