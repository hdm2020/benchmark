#merge experimental and predictive difference value(2 alleles off snp) of TF binding

#read evaluation data
evaldatafile<-'../snpdata/testdata/GVAT_novelbatch_TFSYMBOL.txt'
evaldata<-read.table(evaldatafile,sep='\t',header=T,stringsAsFactors = F)
evaldata<-evaldata[,c('snp','pbs','pval','TF_SYMBOL')]

#JASPAR 2022
modeltf<-'./testdata/evaldata_intermotifbreakR.jaspar587tf.csv'
predfile<-'./testdata/pwm_jaspar2022/results/motifbreakR.jaspar2022.txt'
outfile<-'./testdata/motifbreakR.jaspar2022.merged.expe.pred.results.txt'

#read prediction result
model_tf<-read.csv(modeltf,stringsAsFactors=F)
result<-read.table(predfile,header=T,sep='\t',stringsAsFactors=F)
result<-result[,c('SNP_id','geneSymbol','providerName','alleleDiff','alleleEffectSize')]
colnames(result)[1:3]<-c('snp','TF_SYMBOL','motif')
result$snp<-gsub(':','_',result$snp)
result<-merge(result,model_tf)

#merge and output
df<-merge(evaldata,result,by=c('snp','TF_SYMBOL'))
write.table(df,outfile,col.names=T,row.names=F,quote=F,sep='\t')


#HOCOMOCO v11
modeltf<-'./testdata/evaldata_intermotifbreakR.hocomoco400tf.csv'
predfile<-'./testdata/pwm_hocomocov11/results/motifbreakR.hocomocov11.txt'
outfile<-'./testdata/motifbreakR.hocomocov11.merged.expe.pred.results.txt'

#read prediction result
model_tf<-read.csv(modeltf,stringsAsFactors=F)
result<-read.table(predfile,header=T,sep='\t',stringsAsFactors=F)
result<-result[,c('SNP_id','geneSymbol','providerName','alleleDiff','alleleEffectSize')]
colnames(result)[1:3]<-c('snp','TF_SYMBOL','motif')
result$snp<-gsub(':','_',result$snp)
result[result$TF_SYMBOL=='T','TF_SYMBOL']<-'TBXT'
result<-merge(result,model_tf)

#merge and output
df<-merge(evaldata,result,by=c('snp','TF_SYMBOL'))
write.table(df,outfile,col.names=T,row.names=F,quote=F,sep='\t')
