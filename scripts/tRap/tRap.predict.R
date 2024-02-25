#script for predicting snp's effect on TF binding using tRap

#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="snp file name", metavar="character"),
  make_option(c("-m", "--model"), type="character", default=NULL,
              help="model file name:selected TF motif for predicting", metavar="character"),
  make_option("--pwmdb", type="character", default=NULL,
              help="PWM database", metavar="character"),
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
if (is.null(opt$pwmdb)){
  print_help(opt_parser)
  stop("You need to select a PWM database:jaspar or hocomoco.n", call.=FALSE)
}

#define function
affinity_diff<-function(database,motif_file,seq_wt,seq_mut,out_file){
  ##load package
  if (database=='jaspar'){
    library(tRap,lib.loc = "~/miniconda3/lib/R/library/tRap_jaspar2022/")
    motif_lib<-read.jaspar('./JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt')
    print(paste0('JASPAR data source is loaded:',length(motif_lib)))}
  else if(database=='hocomoco'){
    library(tRap,lib.loc = "~/miniconda3/lib/R/library/tRap_pwmhocomoco11/")
    motif_lib<-read.jaspar('./HOCOMOCOv11_core_HUMAN_mono_jaspar_format.txt')
    print(paste0('HOCOMOCO data source is loaded:',length(motif_lib)))}
  else {
    print('Motif Library is missing')}
  
  library(Biostrings)
  
  #load motif name that need to analyse
  motif<-read.csv(motif_file,stringsAsFactors = F)
  motif_lib<-motif_lib[names(motif_lib) %in% motif$model_name]
  print(paste0('Analysed motif:',length(motif_lib)))
  
  #read SNP sequences(wild and mutation) and subset sequences' name
  SNP_seq_WT <- readDNAStringSet(seq_wt)
  SNP_seq_MUT <- readDNAStringSet(seq_mut)
  
  names(SNP_seq_WT) <- paste(names(SNP_seq_WT),"WT",sep = "_")
  names(SNP_seq_MUT)<- paste(names(SNP_seq_MUT),"MUT",sep = "_")
  
  seq<-c(SNP_seq_WT,SNP_seq_MUT)
  pairs<-data.frame(names(SNP_seq_WT),names(SNP_seq_MUT),stringsAsFactors=F)
  
  #put SNP sequence into a list
  l = list(length=length(seq))
  for (i in 1:length(seq)) {
    l[[i]] = list(desc=names(seq)[i], seq=as.character(seq[i]))
  }
  sequences = l
  print('SNP sequences are loaded')
  
  #modify function:rank.factors.for.pairs
  rank.factors.for.pairs.new <- function(matrices, sequences, pairs){
    desc_names<-c()
    for (i in 1:length(sequences)) {
      desc_names <- append(desc_names,sequences[[i]]$desc)
    }
    dummy <- apply(pairs, 1, function(pair) {
      fa.entry1 = sequences[[which(pair[1] == desc_names)]]
      fa.entry2 = sequences[[which(pair[2] == desc_names)]]
      
      cat(paste("pair:", fa.entry1$desc, fa.entry2$desc, "\n", 
                sep = "\n"))
      pair.res = rank.factors(matrices, fa.entry1$seq, fa.entry2$seq)
      
      results<<-data.frame(seq1 = rep(pair[1], nrow(pair.res)), seq2 = rep(pair[2], nrow(pair.res)), matrix = rownames(pair.res), pair.res)
      results[, "seq1"] = as.character(results[, "seq1"])
      results[, "seq2"] = as.character(results[, "seq2"])
      results[, "matrix"] = as.character(results[, "matrix"])
      write.table(results,out_file,sep="\t", quote=F,row.names = F,append = T)
    })
  }
  
  #calculate affinity differences
  rank.factors.for.pairs.new(motif_lib, sequences, pairs)
}

s1<-paste0(opt$file,'.ref.fa')
s2<-paste0(opt$file,'.alt.fa')
affinity_diff(opt$pwmdb,opt$model,s1,s2,opt$out)
