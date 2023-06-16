# this script can be used to fit gev distribution parameters for new matrices
####modify from ~/handongmei_data/benchmark/tRap/tRap/inst/scripts/fit-gev.R
####use matrix id not matrixindex 

library(tRap)
library(Biostrings)

#####motif PFM matrix####
#data(jaspar)#use hocomoco11 data to replace it
hocomoco11 <- read.jaspar("./HOCOMOCOv11_core_HUMAN_mono_jaspar_format.txt")
print(length(hocomoco11))
print('read hocomoco11 finished')

args = commandArgs()
matrix.path="./motiffiles/"
matrices = read.table(paste0(matrix.path,args[7]),header=F,stringsAsFactors=F)[,1]
cat("matrix index", args[7], matrices, "\n")


####read promoter sequences####
#path of promoter sequences
promoter.path = getwd()
promoter.path = gsub('pwm_hocomocov11','',promoter.path)

# first prepare the sequences
region.size = c(200, 500, 600, 800, 1000, 2000, 5000)

all.seq = list()
for (r in region.size) {
  # load human promoter sequences of the specified length
  #sequences = readFASTA(file.path(promoter.path,#'readFASTA' had been droped in new biostrings package
                                  #paste("Hs.promoters.L", r, ".fa", sep="")))
  #sequences = sapply(sequences, "[[", "seq")
  sequences = readDNAStringSet(file.path(promoter.path,
                                         paste("promoter", r, ".fa", sep="")))
  sequences = sapply(sequences, as.character)
  # replace ? by N
  sequences = gsub("?", "N", sequences, fixed=T)
  sequences = gsub("R", "N", sequences, fixed=T)
  all.seq[[length(all.seq) + 1]] = sequences
}


####fit the gev####
# for all the matrices fit the gev
#load 'fit.gev' function
source('../../tRap/R/fit-gev.R')
params = data.frame()
gev.by.length = data.frame()
for (mat.id in matrices) {
  #mat.id = names(jaspar)[mat.index]
  mat = hocomoco11[[mat.id]]
  print(paste(mat.id,'loaded',sep=' '))
  fit = fit.gev(mat, all.seq)
  params = rbind(params, data.frame(matrix=mat.id, fit$params, stringsAsFactors=F))
  gev.by.length = rbind(gev.by.length, data.frame(matrix=mat.id, fit$gev.by.length, stringsAsFactors=F))
}
#output
output.path="./gev_parameters/"
write.table(params,paste0(output.path,args[7],'.params.txt'),sep = "\t",row.names = F,col.names = T,quote = F)
write.table(gev.by.length,paste0(output.path,args[7],'.length.txt'),sep = "\t",row.names = F,col.names = T,quote = F)
print(paste(args[7],'finished',sep=' '))

