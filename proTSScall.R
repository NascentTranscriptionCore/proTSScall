# usage:
#       Rscript proTSScall.R <tss_list> <read_threshold>

cArgs <- commandArgs()
fileArg <- which(startsWith(cArgs, "--file="))
argsArg <- which(startsWith(cArgs, "--args"))
runningScript <- normalizePath(sapply(strsplit(cArgs[fileArg], "--file="), "[", 2))

tssPath <- cArgs[argsArg+1]
readThreshold <- as.numeric(cArgs[argsArg+2])

tssList <- read.table(tssPath, sep = "\t", row.names = 1, stringsAsFactors = F)

if(startsWith(runningScript, "/n/data1/cores/ntc")){
  scriptsPath <- "/n/data1/cores/ntc/scripts"
} else if(startsWith(runningScript, "/n/data1/hms/bcmp/adelman")){
  scriptsPath <- "/n/data1/hms/bcmp/adelman/Scripts/ntc"
} else{
  print("Warning: invoked script not on 'ntc' or 'adelman' paths")
  scriptsPath <- "/n/data1/cores/ntc/scripts"
}

setwd("bedGraphs")
system("mkdir -p forward reverse")
system("ln -rs *_3pr_forward.bedGraph forward")
system("ln -rs *_3pr_reverse.bedGraph reverse")

outPrefix <- sapply(strsplit(getwd(), "/"), tail, 2)[1]
setwd("forward")
system(paste(paste(scriptsPath, "AdelmanLab/NIH_scripts/bedgraphs2stdBedGraph/bedgraphs2stdBedGraph", sep = "/"), paste(outPrefix, "F", sep = "_")))
setwd("../reverse")
system(paste(paste(scriptsPath, "AdelmanLab/NIH_scripts/bedgraphs2stdBedGraph/bedgraphs2stdBedGraph", sep = "/"), paste(outPrefix, "R", sep = "_")))
setwd("../../")
tssPrefix <- sapply(strsplit(sapply(strsplit(tssPath, "/"), "[", length(strsplit(tssPath, "/")[[1]])), ".tss.txt"), "[", 1)
system("mkdir pro_tss")
system(paste(scriptsPath, "/AdelmanLab/NIH_scripts/make_heatmap/make_heatmap -t 6 -l s -s s --nohead -p bedGraphs/forward/", outPrefix, "_F.bedGraph -m bedGraphs/reverse/", outPrefix, "_R.bedGraph -- ", tssPath, " pro_tss/", outPrefix, "_", tssPrefix, "_25mer_+-2kb.txt -2000 25 160", sep = ""))

tssMat <- read.table(paste("pro_tss/", outPrefix, "_", tssPrefix, "_25mer_+-2kb.txt", sep = ""), header = T, sep = "\t", row.names = 1)
tssReads <- rowSums(tssMat[, which(colnames(tssMat) == "X0.24"):which(colnames(tssMat) == "X125.149")])
tssInactive <- names(tssReads)[tssReads <= readThreshold]
write.table(tssList[tssInactive, ], paste("pro_tss/", outPrefix, "_", tssPrefix, "_inactive_formakeheatmap.txt", sep = ""), quote = F, sep = "\t", col.names = F)
tssActive <- names(tssReads)[tssReads > readThreshold]

# identify 1 dominant TSS per active gene
domTss <- vector("list")

for (i in tssActive){
  if(!(tssList[i, "V7"] %in% names(domTss))){
    domTss[[tssList[i, "V7"]]] <- i
  }
  else{
    # replace with current ENST, if more TSS-proximal reads
    if(tssReads[i] > tssReads[domTss[[tssList[i, "V7"]]]]){
      domTss[[tssList[i, "V7"]]] <- i
    }
    else if (tssReads[i] == tssReads[domTss[[tssList[i, "V7"]]]]){
      if(tssList[i, "V6"] == "+"){
        if(tssList[i, "V4"] < tssList[domTss[[tssList[i, "V7"]]], "V4"]){
          domTss[[tssList[i, "V7"]]] <- i
        }
      }
      else if(tssList[i, "V6"] == "-"){
        if(tssList[i, "V5"] > tssList[domTss[[tssList[i, "V7"]]], "V5"]){
          domTss[[tssList[i, "V7"]]] <- i
        }
      }
      else{
        stop("Unexpected strand symbol")
      }
    }
  }
}

tssNonDom <- tssActive[!(tssActive %in% unlist(domTss))]
write.table(tssList[tssNonDom, ], paste("pro_tss/", outPrefix, "_", tssPrefix, "_non-dominant_formakeheatmap.txt", sep = ""), quote = F, sep = "\t", col.names = F)

uniqueTss <- vector("list")
dupTss <- vector("list")

for (i in unlist(domTss)){
  if(tssList[i, "V6"] == "+"){
    tssStart <- paste(tssList[i, "V3"], tssList[i, "V4"], "+", sep = ".")
  }
  else if(tssList[i, "V6"] == "-"){
    tssStart <- paste(tssList[i, "V3"], tssList[i, "V5"], "-", sep = ".")
  }
  else{
    stop("Unexpected strand symbol")
  }
  
  if(!(tssStart %in% names(uniqueTss))){
    uniqueTss[[tssStart]] <- i
    dupTss[[tssStart]] <- i
  }
  else{
    dupTss[[tssStart]] <- c(dupTss[[tssStart]], i)
    # if TSS exists, take longer transcript
    if(tssList[i, "V10"] > tssList[uniqueTss[[tssStart]], "V10"]){
      uniqueTss[[tssStart]] <- i
    }
    # if TSS exists, and transcript lengths the same, take lower ENSG number
    else if(tssList[i, "V10"] == tssList[uniqueTss[[tssStart]], "V10"]){
      digPos <- regexpr("[[:digit:]]+", c(tssList[i, "V7"], tssList[uniqueTss[[tssStart]], "V7"]))
      geneDigits <- substr(c(tssList[i, "V7"], tssList[uniqueTss[[tssStart]], "V7"]), as.integer(digPos), as.integer(digPos) + (attributes(digPos)$match.length - 1))
      if(which.min(as.numeric(geneDigits)) == 1){
        uniqueTss[[tssStart]] <- i
      }
    }
  }
}

write.table(tssList[rownames(tssList) %in% unlist(uniqueTss), ], paste("pro_tss/", outPrefix, "_", tssPrefix, "_dominant_formakeheatmap.txt", sep = ""), quote = F, sep = "\t", col.names = F)

safPos <- data.frame(chr = character(), startTss = integer(), endTss = integer(), startWhole = integer(), endWhole = integer(), startBody = integer(), endBody = integer(), stringsAsFactors = F)
for (i in names(uniqueTss)){
  loc <- unlist(strsplit(i, "[.]"))
  safPos[i, "chr"] <- loc[1]
  if(loc[3] == "+"){
    safPos[i, "startTss"] <- as.integer(loc[2])
    safPos[i, "endTss"] <- as.integer(loc[2]) + 150
    safPos[i, "startWhole"] <- as.integer(loc[2])
    safPos[i, "endWhole"] <- (as.integer(loc[2]) + tssList[uniqueTss[[i]], "V10"]) - 1
    safPos[i, "startBody"] <- as.integer(loc[2]) + 250
    safPos[i, "endBody"] <- as.integer(loc[2]) + 2250
  }
  else if(loc[3] == "-"){
    safPos[i, "startTss"] <- as.integer(loc[2]) - 150
    safPos[i, "endTss"] <- as.integer(loc[2])
    safPos[i, "startWhole"] <- (as.integer(loc[2]) - tssList[uniqueTss[[i]], "V10"]) + 1
    safPos[i, "endWhole"] <- as.integer(loc[2])
    safPos[i, "startBody"] <- as.integer(loc[2]) - 2250
    safPos[i, "endBody"] <- as.integer(loc[2]) - 250
  }
  else{
    stop("Unexpected strand symbol")
  }
}

posOrderTss <- rownames(safPos)[order(safPos$chr, safPos$startTss)]
posOrderWhole <- rownames(safPos)[order(safPos$chr, safPos$startWhole)]
posOrderBody <- rownames(safPos)[order(safPos$chr, safPos$startBody)]

write.table(cbind(unlist(uniqueTss[posOrderTss]), safPos[posOrderTss, "chr"], safPos[posOrderTss, "startTss"], safPos[posOrderTss, "endTss"], sapply(strsplit(posOrderTss, "[.]"), "[", 3), tssList[unlist(uniqueTss[posOrderTss]), "V8"], sapply(dupTss[posOrderTss], paste, collapse = ";")), paste("pro_tss/", outPrefix, "_", tssPrefix, "_tss.saf", sep = ""), quote = F, sep = "\t", col.names = c("GeneID", "Chr", "Start", "End", "Strand", "GeneName", "AllTx"), row.names = F)
write.table(cbind(unlist(uniqueTss[posOrderWhole]), safPos[posOrderWhole, "chr"], safPos[posOrderWhole, "startWhole"], safPos[posOrderWhole, "endWhole"], sapply(strsplit(posOrderWhole, "[.]"), "[", 3), tssList[unlist(uniqueTss[posOrderWhole]), "V8"], sapply(dupTss[posOrderWhole], paste, collapse = ";")), paste("pro_tss/", outPrefix, "_", tssPrefix, "_whole_gene.saf", sep = ""), quote = F, sep = "\t", col.names = c("GeneID", "Chr", "Start", "End", "Strand", "GeneName", "AllTx"), row.names = F)
write.table(cbind(unlist(uniqueTss[posOrderBody]), safPos[posOrderBody, "chr"], safPos[posOrderBody, "startBody"], safPos[posOrderBody, "endBody"], sapply(strsplit(posOrderBody, "[.]"), "[", 3), tssList[unlist(uniqueTss[posOrderBody]), "V8"], sapply(dupTss[posOrderBody], paste, collapse = ";")), paste("pro_tss/", outPrefix, "_", tssPrefix, "_gene_body.saf", sep = ""), quote = F, sep = "\t", col.names = c("GeneID", "Chr", "Start", "End", "Strand", "GeneName", "AllTx"), row.names = F)
