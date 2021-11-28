library(Biostrings)
library(data.table)

args <- c("pavc.txt","HT")
args <- commandArgs()

pavc <- fread(args[1], header = T, na.strings = "-")
pavc <- pavc[which(pavc$Type != "INV" & pavc$Type != "TRANS")]

refgenome <- readDNAStringSet("refgenome")
qrygenome <- readDNAStringSet("qrygenome")

vcf.header <- data.table(V1 = c("##fileformat=VCFv4.2", "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>", paste0("##contig=<ID=", unique(pavc$Ref_Chr),">")))
vcf <- data.table(V1 = pavc$Ref_Chr, V2 = pavc$Ref_Start, V3 = ".", V4 = ".", V5 = ".", V6 = ".", V7 = ".", V8 = ".", V9 = "GT", V10 = "1/1")
for (i in seq(nrow(vcf))){
    print(i)
    if (!is.na(pavc$Qry_Chr[i])){
        vcf$V4[i] <- as.character(subseq(refgenome[pavc$Ref_Chr[i]], pavc$Ref_Start[i], pavc$Ref_End[i]))
        vcf$V5[i] <- as.character(subseq(qrygenome[pavc$Qry_Chr[i]], pavc$Qry_Start[i], pavc$Qry_End[i]))
    }else{
        vcf$V4[i] <- as.character(subseq(refgenome[pavc$Ref_Chr[i]], pavc$Ref_Start[i], pavc$Ref_End[i]))
        vcf$V5[i] <- as.character(subseq(refgenome[pavc$Ref_Chr[i]], pavc$Ref_Start[i], pavc$Ref_Start[i]))
    }

}
colnames(vcf) <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",args[2])
vcf$ID <- pavc$New_ID
vcf <- vcf[order(vcf$`#CHROM`, vcf$POS)]
fwrite(vcf.header, file = "pavc.vcf", quote = F, row.names = F, col.names = F, sep = "\t")
fwrite(vcf, file = "pavc.vcf", quote = F, row.names = F, col.names = T, sep = "\t", append = T)
