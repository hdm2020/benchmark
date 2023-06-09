#script for predicting snp's effect on TF binding using atSNP 

#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="snp file name", metavar="character"),
  make_option(c("-m", "--model"), type="character", default=NULL,
              help="model file name:selected TF motif for predicting", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default=NULL,
              help="genome name", metavar="character"),
  make_option("--pwmdb", type="character", default=NULL,
              help="PWM database", metavar="character"),
  make_option(c("-n", "--ncore"), type="integer", default=10,
              help="number of cores [default= %default]", metavar="number"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
if (is.null(opt$model)){
  print_help(opt_parser)
  stop("You need to provide some TFs you want to predict.n", call.=FALSE)
}
if (is.null(opt$genome)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input genome).n", call.=FALSE)
}
if (is.null(opt$pwmdb)){
  print_help(opt_parser)
  stop("You need to select a PWM database:jaspar or hocomoco.n", call.=FALSE)
}

#load motif library
library(atSNP)
selectpwm<-function(pwmdb){
  if(pwmdb=='jaspar'){
    pwms<- LoadMotifLibrary(
    filename="/home/handongmei/handongmei_data/benchmark/atSNP/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
    tag = ">", skiprows = 1, skipcols = 1, transpose = TRUE,
    field = 1, sep = c("\t", " ", "\\[", "\\]", ">"),
    pseudocount = 1)
    return(pwms)}
  else if(pwmdb=='hocomoco'){
    pwms<- LoadMotifLibrary(
    filename="/home/handongmei/handongmei_data/benchmark/atSNP/HOCOMOCOv11_core_HUMAN_mono_jaspar_format.txt",
    tag = ">", skiprows = 1, skipcols = 0, transpose = TRUE, field = 1,
    sep = c(">", "\t", " "), pseudocount = 1)
    return(pwms)}
  else {
    print('Motif Library is missing')}
}

pwms<-selectpwm(opt$pwmdb)
tf<-read.csv(opt$model,stringsAsFactors=F)
motif_library<-pwms[tf$model_name]
print(motif_library[1:10])
print('motif library have been loaded')

#load snp data
snpinfo <- LoadSNPData(opt$file, genome.lib = opt$genome, half.window.size = 30, default.par = FALSE, mutation = TRUE)
print('snp data have been loaded')
#default.par = FALSE if number of snps >1000

#Compute affinity scores
atsnp.scores <- ComputeMotifScore(motif_library, snpinfo, ncores = opt$ncore)#ncores:An integer for the number of parallel process. Default: 1
print('scores have been calculated')
atsnp.result <- ComputePValues(motif.lib = motif_library, snp.info = snpinfo,
                               motif.scores = atsnp.scores$motif.scores, ncores = opt$ncore, testing.mc=TRUE)
print('p values have been calculated')
write.table(atsnp.result,opt$out,row.names=F,col.names=T,sep='\t',quote=F)

