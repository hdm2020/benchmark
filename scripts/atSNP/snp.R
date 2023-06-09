#modify snps' format
#five columns:snpid chr snp(position) a1(ref allele) a2(alt allele)
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
colnames(df)<-'snpid'
s1<-function(xx){return(strsplit(xx,'_')[[1]][1])}
s2<-function(xx){return(strsplit(xx,'_')[[1]][2])}
s3<-function(xx){return(strsplit(xx,'_')[[1]][3])}
s4<-function(xx){return(strsplit(xx,'_')[[1]][4])}
df$chr<-apply(df,1,s1)
df$snp<-apply(df,1,s2)
df$a1<-apply(df,1,s3)
df$a2<-apply(df,1,s4)
write.table(df,opt$out,row.names=F,col.names=T,quote=F,sep=' ')

