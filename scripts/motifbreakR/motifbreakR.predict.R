#script for predicting snp's effect on TF binding using atSNP
#Note: for saving time,you can split original snp file into some small files

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
  make_option("--method", type="character", default='default',
              help="method:default or ic or log [default= %default]", metavar="character"),
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
library(motifbreakR)
library(MotifDb)
library(BiocParallel)

selectpwm<-function(pwmdb){
  if(pwmdb=='jaspar'){
    motiflist<-read.csv(opt$model,stringsAsFactors=F)
    xx<-subset(MotifDb,dataSource=='jaspar2022' & organism=='Hsapiens')#691
    motiflibrary<-query(xx,andStrings=c('Hsapiens','jaspar2022'),orStrings=motiflist$model_name)
    return(motiflibrary)
    }
  else if(pwmdb=='hocomoco'){
    motiflist<-read.csv(opt$model,stringsAsFactors=F)
    xx<-subset(MotifDb,dataSource=='HOCOMOCOv11-core-A'|dataSource=='HOCOMOCOv11-core-B'|dataSource=='HOCOMOCOv11-core-C' & organism=='Hsapiens')#400
    motiflibrary<-query(xx,andStrings=c('Hsapiens'),orStrings=motiflist$model_name)
    return(motiflibrary)
    }
  else {
    print('Motif Library is missing')}
}

pwms<-selectpwm(opt$pwmdb)
print('motif library is loaded')
print(length(pwms))

#define function for prediction
motifbreakr_f<-function(inputfile,genome,method,outputfile){
  start<-Sys.time()
  print(start)
  #load snp data from bed file --no rsid
  if(genome=='BSgenome.Hsapiens.UCSC.hg19'){
    library(BSgenome.Hsapiens.UCSC.hg19)
    snps.mb.frombed <- snps.from.file(file = inputfile,
                                      search.genome = BSgenome.Hsapiens.UCSC.hg19,
                                      format = "bed")
  }
  else if(genome=='BSgenome.Hsapiens.UCSC.hg38'){
    library(BSgenome.Hsapiens.UCSC.hg38)
    snps.mb.frombed <- snps.from.file(file = inputfile,
                                      search.genome = BSgenome.Hsapiens.UCSC.hg38,
                                      format = "bed")
  }
  else{
    print('genome is missing')
  }

  print('snp is loaded')
  print(length(snps.mb.frombed))
  
  #run motifbreakR()
  result<-motifbreakR(
    snpList = snps.mb.frombed,
    pwmList = pwms,
    filterp = TRUE,
    threshold = 1,
    method = "default",
    show.neutral = TRUE,
    bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
    BPPARAM = bpparam()
  )
  print('run motifbreakR is finished')
  
  #save result
  df<-data.frame(result)
  df<-df[-10]
  write.table(df,outputfile,row.names=F,col.names=T,quote=F,sep='\t')
  print(paste0(outputfile,' is finished'))
  end<-Sys.time()
  print(end)
  print(end-start) 
}
#run
motifbreakr_f(opt$file,opt$genome,opt$method,opt$out) 

