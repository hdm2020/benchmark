#select every best tf model's merged experimental and predictive difference value of TF binding

#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option("--f1", type="character", default=NULL,
              help="auroc file name", metavar="character"),
  make_option("--f2", type="character", default=NULL,
              help="the file including merged experimental and predictive difference value of TF binding", metavar="character"),
  make_option("--out", type="character", default=NULL,
              help="output file name:merged experimental and predictive difference value of TF binding", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#TFs' auroc
df<-read.table(opt$f1,header=T,sep='\t',stringsAsFactors=F)
#merged experimental and predictive difference value of TF binding
result<-read.table(opt$f2,header=T,sep='\t',stringsAsFactors=F)

result<-result[result$model_name %in% df$model_name,]

write.table(result,opt$out,row.names=F,col.names=T,sep='\t',quote=F)

