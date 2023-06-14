#merge experimental and predictive difference value(2 alleles off snp) of TF binding

#read evaluation data
evaldatafile<-'../snpdata/testdata/GVAT_novelbatch_TFSYMBOL.txt'
evaldata<-read.table(evaldatafile,sep='\t',header=T,stringsAsFactors = F)
evaldata<-evaldata[,c('snp','pbs','pval','TF_SYMBOL')]

#JASPAR 2022
modeltf<-'./testdata/evaldata_interatSNP.jaspar587tf.csv'
predfile<-'./testdata/pwm_jaspar2022/results/atSNP.jaspar2022.txt'
outfile<-'./testdata/atSNP.jaspar2022.merged.expe.pred.results.txt'

#read prediction result
model_tf<-read.csv(modeltf,stringsAsFactors=F)
result<-read.table(predfile,header=T,sep='\t',stringsAsFactors=F)
result<-result[,c('motif','snpid','log_lik_ratio','pval_ref','pval_snp','pval_diff','pval_rank')]
colnames(result)[,1:2]<-c('model_name','snp')
result<-merge(result,model_tf)

#merge and output
df<-merge(evaldata,result,by=c('snp','TF_SYMBOL'))
write.table(df,outfile,col.names=T,row.names=F,quote=F,sep='\t')


#HOCOMOCO v11
modeltf<-'./testdata/evaldata_interatSNP.hocomoco401tf.csv'
predfile<-'./testdata/pwm_hocomocov11/results/atSNP.hocomocov11.txt'
outfile<-'./testdata/atSNP.hocomocov11.merged.expe.pred.results.txt'

#read prediction result
model_tf<-read.csv(modeltf,stringsAsFactors=F)
result<-read.table(predfile,header=T,sep='\t',stringsAsFactors=F)
result<-result[,c('motif','snpid','log_lik_ratio','pval_ref','pval_snp','pval_diff','pval_rank')]
colnames(result)[,1:2]<-c('model_name','TF_SYMBOL')
result<-merge(result,model_tf)

#merge and output
df<-merge(evaldata,result,by=c('snp','TF_SYMBOL'))
write.table(df,outfile,col.names=T,row.names=F,quote=F,sep='\t')
