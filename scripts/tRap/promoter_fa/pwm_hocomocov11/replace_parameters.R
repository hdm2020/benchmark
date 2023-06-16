#replace old parameters-exactgevparams.RData  gevparams.RData sysdata.rda
#gevpaprams.RData
gevparams<-read.table('./summary.params.txt',header=T,sep='\t',stringsAsFactors=F)
rownames(gevparams)<-gevparams$matrix
gevparams<-gevparams[,-1]

#exactgevparams.RData
exactgev<-read.table('./summary.length.txt',header=T,sep='\t',stringsAsFactors=F)
library(tRap)
motif<- read.jaspar("./HOCOMOCOv11_core_HUMAN_mono_jaspar_format.txt")
motif_len<-data.frame(matrix=names(motif))
for (i in motif_len$matrix){motif_len[motif_len$matrix==i,'width']<-ncol(motif[[i]])}
exactgev<-merge(exactgev,motif_len)
exactgev<-exactgev[,c('matrix','width','region.size','loc','scale','shape')]
colnames(exactgev)[3:4]<-c('length','location')

exactgev<-split(exactgev, exactgev[,"length"])
exactgevparams<-list()
for (i in c(200,500,600,800,1000,2000,5000)){
  j<-as.character(i)
  rownames(exactgev[[j]]) <- exactgev[[j]]$matrix
  exactgev[[j]] <- exactgev[[j]][,-1]
  exactgevparams[[i]]<-exactgev[[j]]
}

#save
hocomoco<-motif
save(hocomoco,file="./tRap_pwmhocomoco11/data/hocomoco.RData")
save(gevparams, file="./tRap_pwmhocomoco11/data/gevparams.RData")
save(exactgevparams, file="./tRap_pwmhocomoco11/data/exactgevparams.RData")
save(list=c("gevparams", "exactgevparams"), file="./tRap_pwmhocomoco11/R/sysdata.rda")
