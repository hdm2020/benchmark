#modify snps' format: BED format
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
df$snpid<-gsub('_',':',df$snpid)
s1<-function(xx){return(strsplit(xx,':')[[1]][1])}
s2<-function(xx){return(strsplit(xx,':')[[1]][2])}
df$chr<-apply(df,1,s1)
df$snp<-as.numeric(apply(df,1,s2))
df$start<-df$snp-1
df$score<-0;df$strand<-'+'
df<-df[,c('chr','start','snp','snpid','score','strand')]
write.table(df,opt$out,row.names=F,col.names=F,quote=F,sep='\t')

