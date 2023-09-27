#merge experimental and predictive difference value(2 alleles off snp) of TF binding

library(rhdf5)

modeltf<-'./testdata/evaldata_interSei9315tf.csv'
snpid<-'./testdata/testsnppos.vcf'
evaldatafile<-'../snpdata/testdata/GVAT_novelbatch_TFSYMBOL.txt'
predfile<-'./testdata/results/chromatin-profiles-hdf5/testsnppos_diffs.h5'
outfile<-'./testdata/Sei.merged.expe.pred.results.txt'

#read prediction result
model_tf<-read.csv(modeltf,stringsAsFactors=F)
snp_id<-read.table(snpid,sep='\t',header=F,stringsAsFactors = F)
result<-H5Fopen(predfile)
result<-data.frame(t(result$data))
result$snp<-snp_id$V3
allresult<-data.frame()
for (i in 1:nrow(model_tf)){
  print(i)
  df<-result[,c(ncol(result),model_tf[i,'tf_coord'])]
  colnames(df)<-c('snp','delta_alt_ref')
  df$TF_SYMBOL<-model_tf[i,'TF_SYMBOL']
  df$model_name<-model_tf[i,'model_name']
  allresult<-rbind(allresult,df)
}

rm(result)
#read evaluation data
evaldata<-read.table(evaldatafile,sep='\t',header=T,stringsAsFactors = F)
evaldata<-evaldata[,c('snp','pbs','pval','TF_SYMBOL')]

#merge and output
df<-merge(evaldata,allresult,by=c('snp','TF_SYMBOL'))
write.table(df,outfile,col.names=T,row.names=F,quote=F,sep='\t')

h5closeAll()

