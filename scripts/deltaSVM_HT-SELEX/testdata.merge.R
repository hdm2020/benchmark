#merge experimental and predictive difference value(2 alleles off snp) of TF binding

#read evaluation data
evaldatafile<-'../snpdata/testdata/GVAT_novelbatch_TFSYMBOL.txt'
evaldata<-read.table(evaldatafile,sep='\t',header=T,stringsAsFactors = F)
evaldata<-evaldata[,c('snp','pbs','pval','TF_SYMBOL')]

outfile<-'./testdata/deltaSVM_HT-SELEX.merged.expe.pred.results.txt'

#94 high confidence TFs
modeltf<-'./testdata/evaldata_inter94tf.csv'
predfile<-'./testdata/94tf_results/testdata.pbs.pred.tsv'

#read prediction result
model_tf<-read.csv(modeltf,stringsAsFactors=F)
result<-read.table(predfile,header=T,stringsAsFactors=F,sep='\t')[,c('snp','deltaSVM','TF')]
result[result$TF=='T','TF']<-'TBXT'
colnames(result)[3]<-'TF_SYMBOL'
result<-merge(result,model_tf[,c('TF_SYMBOL','model_name')])
result<-result[,c('snp','deltaSVM','TF_SYMBOL','model_name')]

df1<-merge(evaldata,result,by=c('snp','TF_SYMBOL'))
#write.table(df1,'./testdata/deltaSVM_HT-SELEX.94tf.merged.expe.pred.results.txt',col.names=T,row.names=F,quote=F,sep='\t')

#other 439 TFs
modeltf<-'./testdata/evaldata_inter439tf.csv'

#read prediction result
model_tf<-read.csv(modeltf,stringsAsFactors=F)
result<-data.frame()
for (i in 1:nrow(model_tf)){
  df<-read.table(paste0('testdata/439tf_results/',model_tf[i,'model_name'],'.deltasvm_scores.tsv'),header=F,stringsAsFactors=F,sep='\t')
  colnames(df)<-c('snp','deltaSVM')
  df$TF_SYMBOL<-model_tf[i,'TF_SYMBOL']
  df$model_name<-model_tf[i,'model_name']
  result<-rbind(result,df)
}

df2<-merge(evaldata,result,by=c('snp','TF_SYMBOL'))
#write.table(df2,'./testdata/deltaSVM_HT-SELEX.439tf.merged.expe.pred.results.txt',col.names=T,row.names=F,quote=F,sep='\t')

#merge and output
alldf<-rbind(df1,df2)
write.table(alldf,outfile,col.names=T,row.names=F,quote=F,sep='\t')

