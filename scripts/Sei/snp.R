#modify snps' format
#five columns: chr pos id ref alt
#!/usr/bin/env Rscript
library("optparse")
 
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="snp file name", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

#script
df<-read.table(opt$file,header=F,stringsAsFactors=F)
colnames(df)<-'id'
s1<-function(xx){return(strsplit(xx,'_')[[1]][1])}
s2<-function(xx){return(strsplit(xx,'_')[[1]][2])}
s3<-function(xx){return(strsplit(xx,'_')[[1]][3])}
s4<-function(xx){return(strsplit(xx,'_')[[1]][4])}
df$chr<-apply(df,1,s1)
df$pos<-apply(df,1,s2)
df$ref<-apply(df,1,s3)
df$alt<-apply(df,1,s4)
write.table(df[,c('chr','pos','id','ref','alt')],opt$out,row.names=F,col.names=F,quote=F,sep='\t')

