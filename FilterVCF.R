library(data.table)

args <- c("pavc.norm.vcf","50")
args <- commandArgs(T)

len <- as.integer(args[2])
vcf <- fread(args[1], header = T, sep = "\t")
vcf$MARK <- "keep"
vcf$reflen <- nchar(vcf$REF)
vcf$altlen <- nchar(vcf$ALT)
vcf$MARK[which(vcf$reflen < len & vcf$altlen < len)] <- "remove"
vcf.keep <- vcf[which(vcf$MARK == "keep")]
vcf.keep[,c(11,12,13)] <- NULL

vcf.header <- data.table(V1 = c("##fileformat=VCFv4.2", "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>", paste0("##contig=<ID=", unique(pavc$Ref_Chr),">")))
fwrite(vcf.header, file = paste0("pavc.norm.FilterLen",args[2],".vcf"), quote = F, row.names = F, col.names = F, sep = "\t")
fwrite(vcf.keep, file = paste0("pavc.norm.FilterLen",args[2],".vcf"), quote = F, row.names = F, col.names = T, sep = "\t", append = T)
