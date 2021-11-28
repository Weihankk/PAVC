library(data.table)

#args <- c("syri.out")
args <- commandArgs()

syri.out <- fread(args[1], header = F)

# Step 1. Extract DSTRs and sort them according to the refgenome and qrygenome, respectively.
dstr <- syri.out[which(syri.out$V11 %in% c("SYN", "INV", "TRANS", "INVTR", "DUP", "INVDP", "NOTAL"))]
dstr$V2 <- as.integer(dstr$V2)
dstr$V3 <- as.integer(dstr$V3)
dstr$V7 <- as.integer(dstr$V7)
dstr$V8 <- as.integer(dstr$V8)
dstr.onref <- dstr[which(dstr$V1 != "-" & dstr$V12 != "copygain")]
dstr.onqry <- dstr[which(dstr$V6 != "-" & dstr$V12 != "copyloss")]
dstr.onref <- dstr.onref[order(dstr.onref$V1, dstr.onref$V2, dstr.onref$V3)]
dstr.onqry <- dstr.onqry[order(dstr.onqry$V6, dstr.onqry$V7, dstr.onqry$V8)]

# Step 2. Classify some DSTRs into presence/absence.
dstr.onqry$Classify <- NA
dstr.onref$Classify <- NA
dstr.onref$Classify[which(dstr.onref$V11 != "SYN" & dstr.onref$V11 != "INV" & dstr.onref$V11 != "TRANS" & dstr.onref$V11 != "INVTR")] <- "Absence"
dstr.onqry$Classify[which(dstr.onqry$V11 != "SYN" & dstr.onqry$V11 != "INV" & dstr.onqry$V11 != "TRANS" & dstr.onqry$V11 != "INVTR")] <- "Presence"

# Step 3. Detect the breakpoints of these PAVs on refgenome.
ref.chr <- unique(dstr.onref$V1)
ref.list <- list()
for (chr in ref.chr){
    dstr.onref.chr <- dstr.onref[which(dstr.onref$V1 == chr)]
    merge.chr <- chr
    merge.start <- c()
    merge.end <- c()
    merge.id <- c()
    merge.mark <- "Syntenic"
    
    
    for (i in seq(nrow(dstr.onref.chr))){
        if (is.na(dstr.onref.chr$Classify[i])){
            
            if (merge.mark == "Absence"){
                # Complete merge a block
                tmp.dt <- data.table(V1 = merge.chr, V2 = min(merge.start), V3 = max(merge.end), V4 = "-", V5 = "-", V6 = NA, V7 = NA, V8 = NA, V9 = paste(merge.id, collapse = "|"), V10 = "-", V11 = "MERGE", V12 = "-", Classify = "Absence")
                ref.list[[paste(chr,i-1)]] <- tmp.dt
                ref.list[[paste(chr,i)]] <- dstr.onref.chr[i]
                merge.mark <- "Syntenic"
                merge.start <- c()
                merge.end <- c()
                merge.id <- c()
                
            }else{
                ref.list[[paste(chr,i)]] <- dstr.onref.chr[i]
            }
    
        }else{
            merge.start <- c(merge.start, dstr.onref.chr$V2[i])
            merge.end <- c(merge.start, dstr.onref.chr$V3[i])
            merge.id <- c(merge.id, dstr.onref.chr$V9[i])
            merge.mark <- "Absence"
        }
    }
}
dstr.onref.merge <- rbindlist(ref.list)

qry.chr <- unique(dstr.onqry$V6)
qry.list <- list()
for (chr in qry.chr){
    dstr.onqry.chr <- dstr.onqry[which(dstr.onqry$V6 == chr)]
    merge.chr <- chr
    merge.start <- c()
    merge.end <- c()
    merge.id <- c()
    merge.mark <- "Syntenic"
    
    
    for (i in seq(nrow(dstr.onqry.chr))){
        if (is.na(dstr.onqry.chr$Classify[i])){
            
            if (merge.mark == "Absence"){
                # Complete merge a block
                tmp.dt <- data.table(V1 = dstr.onqry.chr$V1[i], V2 = dstr.onqry.chr$V2[i], V3 = dstr.onqry.chr$V2[i], V4 = "-", V5 = "-", V6 = merge.chr, V7 = min(merge.start), V8 = max(merge.end), V9 = paste(merge.id, collapse = "|"), V10 = "-", V11 = "MERGE", V12 = "-", Classify = "Presence")
                qry.list[[paste(chr,i-1)]] <- tmp.dt
                qry.list[[paste(chr,i)]] <- dstr.onqry.chr[i]
                merge.mark <- "Syntenic"
                merge.start <- c()
                merge.end <- c()
                merge.id <- c()
                
            }else{
                ref.list[[paste(chr,i)]] <- dstr.onref.chr[i]
            }
            
        }else{
            merge.start <- c(merge.start, dstr.onqry.chr$V7[i])
            merge.end <- c(merge.start, dstr.onqry.chr$V8[i])
            merge.id <- c(merge.id, dstr.onqry.chr$V9[i])
            merge.mark <- "Absence"
        }
    }
}
dstr.onqry.merge <- rbindlist(qry.list)

dstr.onref.onqry <- rbind(dstr.onref.merge, dstr.onqry.merge[which(dstr.onqry.merge$V11 != "SYN")])
dstr.onref.onqry.sort <- dstr.onref.onqry[order(dstr.onref.onqry$V1, dstr.onref.onqry$V2, dstr.onref.onqry$V3)]

## Step 4. Classify the DSEQs into presence/absence.
dstr.syn <- dstr.onref.onqry[which(dstr.onref.onqry$V11 == "SYN")]
syn.list <- list()
for (i in seq(nrow(dstr.syn))){
    syn.name <- dstr.syn$V9[i]
    syn.dseq <- syri.out[which(syri.out$V10 == syn.name)]
    syn.dseq$V2 <- as.integer(syn.dseq$V2)
    syn.dseq$V3 <- as.integer(syn.dseq$V3)
    syn.dseq$V7 <- as.integer(syn.dseq$V7)
    syn.dseq$V8 <- as.integer(syn.dseq$V8)
    if (nrow(syn.dseq) > 0){
        for (s in seq(nrow(syn.dseq))){
            dseq.type <- syn.dseq$V11[s]
            if (dseq.type == "DEL"){
                # DEL -- Absence
                tmp.dt <- data.table(V1 = syn.dseq$V1[s], V2 = syn.dseq$V2[s], V3 = syn.dseq$V3[s], V4 = "-", V5 = "-", V6 = syn.dseq$V6[s], V7 = syn.dseq$V7[s], V8 = syn.dseq$V8[s], V9 = syn.dseq$V9[s], V10 = syn.dseq$V10[s], V11 = syn.dseq$V11[s], V12 = "-", Classify = "Absence")
                syn.list[[syn.dseq$V9[s]]] <- tmp.dt
            }else if (dseq.type == "INS"){
                # INS -- Presence
                tmp.dt <- data.table(V1 = syn.dseq$V1[s], V2 = syn.dseq$V2[s], V3 = syn.dseq$V3[s], V4 = "-", V5 = "-", V6 = syn.dseq$V6[s], V7 = syn.dseq$V7[s], V8 = syn.dseq$V8[s], V9 = syn.dseq$V9[s], V10 = syn.dseq$V10[s], V11 = syn.dseq$V11[s], V12 ="-", Classify = "Presence")
                syn.list[[syn.dseq$V9[s]]] <- tmp.dt
            }else if (dseq.type == "CPG"){
                # CPG -- Presence
                tmp.dt <- data.table(V1 = syn.dseq$V1[s], V2 = syn.dseq$V2[s], V3 = syn.dseq$V3[s], V4 = "-", V5 = "-", V6 = syn.dseq$V6[s], V7 = syn.dseq$V7[s] + 1, V8 = syn.dseq$V8[s] - 1, V9 = syn.dseq$V9[s], V10 = syn.dseq$V10[s], V11 = syn.dseq$V11[s], V12 ="-", Classify = "Presence")
                syn.list[[syn.dseq$V9[s]]] <- tmp.dt
            }else if (dseq.type == "CPL"){
                # CPL -- Absence
                tmp.dt <- data.table(V1 = syn.dseq$V1[s], V2 = syn.dseq$V2[s] + 1, V3 = syn.dseq$V3[s] - 1, V4 = "-", V5 = "-", V6 = syn.dseq$V6[s], V7 = syn.dseq$V7[s], V8 = syn.dseq$V8[s], V9 = syn.dseq$V9[s], V10 = syn.dseq$V10[s], V11 = syn.dseq$V11[s], V12 ="-", Classify = "Absence")
                syn.list[[syn.dseq$V9[s]]] <- tmp.dt
            }else if (dseq.type == "HDR"){
                # HDR in refgenome is absence
                # HDR in qrygenome is presence
                hdr.absence <- data.table(V1 = syn.dseq$V1[s], V2 = syn.dseq$V2[s] - 1, V3 = syn.dseq$V3[s], V4 = "-", V5 = "-", V6 = syn.dseq$V6[s], V7 = syn.dseq$V7[s] - 1, V8 = syn.dseq$V7[s] - 1, V9 = syn.dseq$V9[s], V10 = syn.dseq$V10[s], V11 = syn.dseq$V11[s], V12 ="-", Classify = "Absence")
                hdr.presence <- data.table(V1 = syn.dseq$V1[s], V2 = syn.dseq$V3[s] + 1, V3 = syn.dseq$V3[s] + 1, V4 = "-", V5 = "-", V6 = syn.dseq$V6[s], V7 = syn.dseq$V7[s], V8 = syn.dseq$V8[s] + 1, V9 = syn.dseq$V9[s], V10 = syn.dseq$V10[s], V11 = syn.dseq$V11[s], V12 ="-", Classify = "Presence")
                tmp.dt <- rbind(hdr.absence, hdr.presence)
                syn.list[[syn.dseq$V9[s]]] <- tmp.dt
            }else if (dseq.type == "TDM"){
                tmp.dt <- data.table(V1 = syn.dseq$V1[s], V2 = syn.dseq$V2[s], V3 = syn.dseq$V3[s], V4 = "-", V5 = "-", V6 = syn.dseq$V6[s], V7 = syn.dseq$V7[s], V8 = syn.dseq$V8[s], V9 = syn.dseq$V9[s], V10 = syn.dseq$V10[s], V11 = syn.dseq$V11[s], V12 = "-", Classify = "NA")
                if ((tmp.dt$V3 - tmp.dt$V2) > (tmp.dt$V8 - tmp.dt$V7)){
                    # This is a absence
                    tmp.dt$Classify <- "Absence"
                }else{
                    # This is a presence
                    tmp.dt$Classify <- "Presence"
                }
                syn.list[[syn.dseq$V9[s]]] <- tmp.dt
            }
        }
    }
}

syn.merge <- rbindlist(syn.list)

#### Step 5. Combine PAV
dstr.dseq <- list()
for (i in seq(nrow(dstr.onref.onqry.sort))){
    if (dstr.onref.onqry.sort$V11[i] == "SYN"){
        dstr.dseq[[dstr.onref.onqry.sort$V9[i]]] <- syn.merge[which(syn.merge$V10 == dstr.onref.onqry.sort$V9[i])]
    }else{
        dstr.dseq[[dstr.onref.onqry.sort$V9[i]]] <- dstr.onref.onqry.sort[i]
    }
}

dstr.dseq <- rbindlist(dstr.dseq)
dstr.dseq$Classify[which(dstr.dseq$V11 == "INV")] <- "INV"
dstr.dseq$Classify[which(dstr.dseq$V11 == "TRANS" | dstr.dseq$V11 == "INVTR")] <- "TRANS"

Final.SyRI <- data.table(Ref_Chr = dstr.dseq$V1, Ref_Start = dstr.dseq$V2, Ref_End = dstr.dseq$V3,
                         Qry_Chr = dstr.dseq$V6, Qry_Start = dstr.dseq$V7, Qry_End = dstr.dseq$V8,
                         New_ID = NA,
                         Type = dstr.dseq$Classify, SyRI_Type = dstr.dseq$V11, SyRI_ID = dstr.dseq$V9)

Final.SyRI <- Final.SyRI[order(Final.SyRI$Ref_Chr, Final.SyRI$Ref_Start)]
Final.SyRI$New_ID <- paste0(Final.SyRI$Type, seq(nrow(Final.SyRI)))

fwrite(Final.SyRI, file = "pavc.txt", quote = F, row.names = F, col.names = T, sep = "\t", na = "-")
