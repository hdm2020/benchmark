#merge experimental and predictive difference value(2 alleles of snp) of TF binding

#read evaluation data
evaldatafile<-'../../snpdata/testdata/GVAT_novelbatch_TFSYMBOL.txt'
evaldata<-read.table(evaldatafile,sep='\t',header=T,stringsAsFactors = F)
evaldata<-evaldata[,c('snp','pbs','pval','TF_SYMBOL')]

#JASPAR 2022
modeltf<-'./testdata/evaldata_intertRap.jaspar587tf.csv'
predfile<-'./testdata/pwm_jaspar2022/results/tRap.jaspar2022.txt'
outfile<-'./testdata/tRap.jaspar2022.merged.expe.pred.results.txt'

#read prediction result
model_tf<-read.csv(modeltf,stringsAsFactors=F)
result<-read.table(predfile,header=T,sep='\t',stringsAsFactors=F)
result<-result[,c(1,3:11)]
result<-result[result$seq1!='seq1',]
colnames(result)[1:2]<-c('snp','model_name')
result$snp<-gsub('_WT','',result$snp)
result<-merge(result,model_tf)

#merge and output
df<-merge(evaldata,result,by=c('snp','TF_SYMBOL'))
write.table(df,outfile,col.names=T,row.names=F,quote=F,sep='\t')


#HOCOMOCO v11
modeltf<-'./testdata/evaldata_intertRap.hocomoco401tf.csv'
predfile<-'./testdata/pwm_hocomocov11/results/tRap.hocomocov11.txt'
outfile<-'./testdata/tRap.hocomocov11.merged.expe.pred.results.txt'

#read prediction result
model_tf<-read.csv(modeltf,stringsAsFactors=F)
result<-read.table(predfile,header=T,sep='\t',stringsAsFactors=F)
result<-result[,c(1,3:11)]
result<-result[result$seq1!='seq1',]
colnames(result)[1:2]<-c('snp','model_name')
result$snp<-gsub('_WT','',result$snp)
result<-merge(result,model_tf)

#merge and output
df<-merge(evaldata,result,by=c('snp','TF_SYMBOL'))
write.table(df,outfile,col.names=T,row.names=F,quote=F,sep='\t')

