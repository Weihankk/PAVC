library(data.table)

args <- c("pavc.norm.FilterLen50.vcf","pavc.txt")
args <- commandArgs(T)

vcf <- fread(args[1])
pavc <- fread(args[2])
other <- pavc[which(pavc$Type != "Presence" & pavc$Type != "Absence")]


vcf$REFLEN <- nchar(vcf$REF)
vcf$ALTLEN <- nchar(vcf$ALT)
vcf$LEN <- vcf$REFLEN
vcf$LEN[which(vcf$ALTLEN > vcf$REFLEN)] <- vcf$ALTLEN[which(vcf$ALTLEN > vcf$REFLEN)] 

other$Qry_Start <- as.integer(other$Qry_Start)
other$Qry_End <- as.integer(other$Qry_End)
other$REFLEN <- other$Ref_End - other$Ref_Start
other$ALTLEN <- other$Qry_End - other$Qry_Start
other$LEN <- other$REFLEN
other$LEN[which(other$ALTLEN > other$REFLEN)] <- other$ALTLEN[which(other$ALTLEN > other$REFLEN)] 

res <- data.table(V1 = c("Presence count:", "Absence count:", "Translocation count:", "Inversion count:",
                         "Presence length(bp):", "Absence length(bp):", "Translocation length(bp):", "Inversion length(bp):",
                         "Total count:", "Total length(bp):"),
                  V2 = c(nrow(vcf[which(vcf$ALTLEN > vcf$REFLEN)]), nrow(vcf[which(vcf$ALTLEN < vcf$REFLEN)]), nrow(other[which(other$Type == "TRANS")]), nrow(other[which(other$Type == "INV")]),
                         sum(vcf$LEN[which(vcf$ALTLEN > vcf$REFLEN)]), sum(vcf$LEN[which(vcf$ALTLEN < vcf$REFLEN)]), sum(other$LEN[which(other$Type == "TRANS")]), sum(other$LEN[which(other$Type == "INV")]),
                         nrow(vcf) + nrow(other), sum(vcf$LEN) + sum(other$LEN)))
fwrite(res, file = "pavc.final.stat.txt", sep = "\t", quote = F, col.names = F, row.names = F)
