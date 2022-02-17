####-------------------------library packages---------------------------####
library(GSVA)
library(GSEABase)
library(cluster)
library(seqinr)
library(ape)
library(msa)
library(sva)
library(MAESTRO)
library(Seurat)
library(jsonlite)
library(facetscales)
library(data.table)
library(ggplot2)
library(dplyr)
library(plyr)
library(reshape2)
library(abind)
library(stringr)
library(RColorBrewer)
library(ggpubr)
library(scales)
library(ggrepel)
library(tidyr)
library(GGally)
library(KEGG.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(pheatmap)
library(ggfortify)
library(edgeR)
library(ggforce)
library(ggdendro)
library(grid)
library(maftools)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(viridis)
library(ggridges)
library(ggsignif)
library(ggalt)
library(ppcor)
library(ggExtra)


#### read in meta information
patient.info <- read.table(file = '~/Documents/BCR_Project_Xiaoqi/Varscan_mutaion_calling/patient_labels.csv',sep = ",",header = TRUE)
pat.info <- melt(patient.info,id.vars = "patientID")
pat.diagnosis <- sapply(c(1:122), function(x){
  name <- paste("FZ-",x,sep = "")
  if(name %in% as.character(patient.info$normal)){return ("Normal")}
  else if(name %in% as.character(patient.info$tumor)) {return ("Tumor")}
  else{return ("")}
})
pat.diagnosis <- pat.diagnosis[pat.diagnosis != ""]
##all patient information including unpaired samples
all.pat.diagnosis <- sapply(c(1:122), function(x){
  name <- paste("FZ-",x,sep = "")
  if(name %in% as.character(patient.info$normal)){return ("Normal")}
  else {return ("Tumor")}
})
###group paired samples
normal <- c(1:18,41,42:70,112,72:81)
tumor <- c(19:24,27:30,33:40,82,83:111,71,113:122)
t.sample <- paste("FZ.",tumor,sep = "")
n.sample <- paste("FZ.",normal,sep = "")

#### read in TPM data
TPM <- read.table(file = "~/Documents/BCR_Project_Xiaoqi/20171024-RNA-Seq data/lung_TPM_all.csv",sep = ",",header = TRUE,row.names = 1)
t.TPM <- log2(TPM+1)

#### read in immune reperotire processed data
# for BCR
load(file = "~/Documents/Description_study/data/TRUST3_MIA_BCRdata_118.Rdata")
load(file = "~/Documents/Description_study/data/TRUST3_MIA_cluster_processYYC_118.Rdata")
load(file = "~/Documents/Description_study/data/TRUST3_MIA_BCR_clonality_118.Rdata")
load(file = "~/Documents/Description_study/data/TRUST3_MIA_BCR_SHM_118.Rdata")
load(file = "~/Documents/Description_study/data/TRUST3_MIA_BCR_lib_reads.Rdata")
# for TCR
load(file = "~/Documents/Description_study/data/TRUST3_MIA_TCRdata_118.Rdata")
load(file = "~/Documents/Description_study/data/TRUST3_MIA_TCR_lib_reads.Rdata")
load(file = "~/Documents/Description_study/data/TRUST3_MIA_TCR_clonality_118.Rdata")
# Ig isotypes
Igs <- c('IGHD','IGHM','IGHE','IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4')
Igs.order <- c("IGHM","IGHD","IGHG3","IGHG1","IGHA1","IGHG2","IGHG4","IGHE","IGHA2")


####function of GO enrichment
GSEAGO <- function(geneList,ont,minGSSize,nPerm,pcut,tsize,msize,title){
  gsea.go <- gseGO(geneList = geneList, 
                   OrgDb = org.Hs.eg.db, 
                   keytype = "ENTREZID",
                   ont = ont, 
                   nPerm = nPerm, # number permutations
                   minGSSize = minGSSize,
                   pvalueCutoff = 1, # padj cutoff value
                   verbose = FALSE)
  kegg <- gsea.go@result %>%
    #mutate(Count = paste("n = ",str_count(core_enrichment,"/")+1, sep = "")) %>%
    mutate(group = ifelse(NES > 0,"Up-regulated","Down-regulated")) %>%
    mutate(significance = cut(p.adjust,breaks=c(-Inf, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ""))) %>%
    dplyr::filter(p.adjust < pcut) %>%
    arrange(desc(NES))
  kegg.res <- NULL
  if(dim(subset(kegg, NES > 0))[1] >= 10){
    kegg.res <- rbind(kegg.res,kegg[1:10,])
  }
  if(dim(subset(kegg, NES > 0))[1] < 10){
    kegg.res <- rbind(kegg.res,subset(kegg, NES > 0))
  }
  if(dim(subset(kegg, NES < 0))[1] >= 10){
    total <- dim(kegg)[1]
    kegg.res <- rbind(kegg.res,kegg[(total-9):total,])
  }
  if(dim(subset(kegg, NES < 0))[1] < 10){
    kegg.res <- rbind(kegg.res,subset(kegg, NES < 0))
  }
  kegg.res$Description <- factor(kegg.res$Description,levels = levels(reorder(kegg.res$Description,kegg.res$NES)))
  kegg.res.new <- kegg.res[match(levels(kegg.res$Description),kegg.res$Description),]
  return(kegg.res.new)
}

####function of KEGG pathway enrichment
GSEAKEGG <- function(geneList,minGSSize,nPerm,pcut,tsize,msize,title,checkgenes){
  GeneIDSymbol <- toTable(org.Hs.egSYMBOL)
  gsea.kegg <- gseKEGG(geneList = geneList, 
                       organism = "hsa", 
                       keyType = "kegg", 
                       nPerm = nPerm, # number permutations
                       minGSSize = minGSSize,
                       pvalueCutoff = 1, # padj cutoff value
                       verbose = FALSE)
  kegg.res <- gsea.kegg@result %>%
    #mutate(Count = paste("n = ",str_count(core_enrichment,"/")+1, sep = "")) %>%
    # mutate(Count = sapply(core_enrichment, function(x)
    #   paste("n = ",length(intersect(checkgenes,subset(GeneIDSymbol,gene_id %in% unlist(strsplit(x, split="/")), select = "symbol")[,1])),sep = ""))) %>%
    mutate(group = ifelse(NES > 0,"Up-regulated","Down-regulated")) %>%
    mutate(significance = cut(p.adjust,breaks=c(-Inf, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ""))) %>%
    dplyr::filter(p.adjust <= pcut) 
  kegg.res$Description <- factor(kegg.res$Description,levels = levels(reorder(kegg.res$Description,kegg.res$NES)))
  kegg.res.new <- kegg.res[match(levels(kegg.res$Description),kegg.res$Description),]
  return(kegg.res.new)
}

####function of barplot for GSEA
GSEAPlot <- function(kegg.res.new,category){
  max.nes <- ceiling(max(kegg.res.new$NES[which(kegg.res.new$NES > 0)]))
  min.nes <- -ceiling(abs(min(kegg.res.new$NES[which(kegg.res.new$NES < 0)])))
  nes.pos <- length(which(kegg.res.new$NES > 0))
  nes.neg <- length(which(kegg.res.new$NES < 0))
  ggplot(kegg.res.new,aes(x=Description,y=NES,fill=as.factor(group))) + 
    xlab(category)+
    ggtitle(title)+
    geom_bar(stat='identity') +
    theme_bw()+
    theme(#axis.text.x=element_blank(),
      plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x=element_text(size=12,face = "bold",hjust=1),
      axis.title.x = element_text(size = 12,face = "bold"),
      axis.title.y = element_text(size = 12,face = "bold"),
      legend.position='bottom',
      panel.background = element_rect(colour = "black", size=1))+
    scale_fill_manual(name = "", values = c("Up-regulated" = "#E41A1C", "Down-regulated" = "#377EB8")) +
    annotate("text",label=kegg.res.new$Description[which(kegg.res.new$NES < 0)],
             x=c(1:nes.neg),
             y=rep(0,nes.neg),size=tsize,hjust=0)+  #negative score
    annotate("text",label=kegg.res.new$Description[which(kegg.res.new$NES > 0)],
             x=c((nes.neg+1):dim(kegg.res.new)[1]),
             y=rep(0,nes.pos),size=tsize,hjust=1)+  #positive score
    # annotate("text",label=kegg.res.new$Count,
    #          x=c(1:dim(kegg.res.new)[1]),
    #          y=c(rep(min.nes/2,nes.neg),rep(max.nes/2,nes.pos)),size=tsize)+
    annotate("text",label=kegg.res.new$significance,
             x=c(1:dim(kegg.res.new)[1]),
             y=c(rep(min.nes,nes.neg),rep(max.nes,nes.pos)),size=msize,color = "black")+
    coord_flip(expand = TRUE)
}

##function of comparing signature gene expression between tumor and normal
GeneCompare <- function(gene,TPM.dat,mark){
  g <- TPM.dat[,c(gene,"Tissue")]
  colnames(g) <- c("expression","Tissue")
  p <- signif(t.test(expression ~ Tissue,data = g)$p.value,4)
  text <- as.character(cut(p, right = FALSE, breaks = c(0, 0.001, 0.01, 0.05,1), labels = c("***","**","*"," ")))
  if(mark == "paired"){
    g.dat <- g %>% arrange(Tissue) %>% mutate(Group = rep(c(1:(dim(g)[1]/2)),times = 2))
    g.dat$Tissue <- factor(g.dat$Tissue, levels = c("Tumor","Normal"))
    gplot <- ggplot(g.dat, aes(Tissue, expression,fill = Tissue)) +
      geom_boxplot() + 
      geom_point(size = 0.5)+
      geom_line(aes(x = Tissue, y = expression, group = Group),color = "grey")+
      theme_bw() +
      #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      scale_fill_manual(values = c("#E41A1C","#377EB8"))+
      ggtitle(gene) +
      xlab("")+
      ylab("mRNA expression")+
      #annotate("text",label = paste("P = ",p,sep = ''),x = 1.5, y = max(g$expression))+
      annotate("text",label = text,x = 1.5, y = max(g$expression)*0.95,size = 7, color = "red")+
      theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
            axis.text.x=element_text(angle=0,size=12,face = "bold",margin = margin(t = 0, b = 10),hjust = 0.5),
            axis.text.y = element_text(size = 12,face = "bold"),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold"),
            legend.position = "none") 
    print (gene)
  }
  if(mark == "unpaired"){
    gplot <- ggplot(g, aes(Tissue, expression,fill = Tissue)) +
      geom_boxplot() + 
      geom_point(size = 0.5)+
      #geom_line(aes(x = Tissue, y = expression, group = Group),color = "grey")+
      theme_bw() +
      #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      scale_fill_manual(values = c("#377EB8","#E41A1C"))+
      ggtitle(gene) +
      xlab("")+
      ylab("mRNA expression")+
      #annotate("text",label = paste("P = ",p,sep = ''),x = 1.5, y = max(g$expression))+
      annotate("text",label = text,x = 1.5, y = max(g$expression)*0.95,size = 6.5, color = "red")+
      theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
            axis.text.x=element_text(angle=0,size=12,face = "bold",margin = margin(t = 0, b = 10),hjust = 0.5),
            axis.text.y = element_text(size = 2,face = "bold"),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold"),
            legend.position = "none") 
    print (gene)
  }
  return (gplot)
}

##function of comparing infiltration in tumor and normal
CompareInfiltration <- function(Infil.est,n,pat.diagnosis,pair,title,tool,xsize){
  infil.melt <- melt(Infil.est,id.vars = "sampleID")
  colnames(infil.melt) <- c("ID","Immune.cell","Estimation")
  infil.melt$Group <- c(rep(pat.diagnosis,n))
  ##calculate significance 
  Infil.dat <- Infil.est
  pat.new <- data.frame(row.names = Infil.est$sampleID, pat.diagnosis)
  Infil.dat$Group <- pat.new$pat.diagnosis
  all.p <- NULL
  for(i in colnames(Infil.dat)[2:(n+1)]){
    subdat <- subset(Infil.dat,select = c(i,"Group"))
    colnames(subdat) <- c("cell","group")
    P <- wilcox.test(cell~group,data=subdat,paired = pair)$p.value
    all.p <- c(all.p,P)
  }
  p.star <- cut(all.p,breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "ns"))
  ##ggplot
  ggplot(infil.melt,aes(x=Immune.cell, y=Estimation,fill=Group))+
    geom_boxplot(outlier.size = 0.5)+
    ggtitle(title)+ 
    ylab("Immune cell Abundance")+
    scale_fill_manual(values = c("#377EB8","#E41A1C"))+
    annotate("text",x = 1:n, y = as.numeric(apply(Infil.est[2:(n+1)],2,max))+0.1, label = p.star, color = "red",size = 7)+
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
          axis.text.x=element_text(angle=30,size=xsize,face = "bold",hjust=1),
          axis.text.y = element_text(size = 12,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12,face = "bold"),
          legend.text = element_text(size = 12),
          legend.position='bottom')
  
}

##function of comparing infiltration in different stages
CompareInfiltrationStage <- function(Infil.est,n,pat.diagnosis,pair,title,xsize){
  infil.melt <- melt(Infil.est,id.vars = "sampleID")
  colnames(infil.melt) <- c("ID","Immune.cell","Estimation")
  infil.melt$Group <- c(rep(pat.diagnosis,n))
  ##calculate significance 
  Infil.dat <- Infil.est
  pat.new <- data.frame(row.names = Infil.est$sampleID, pat.diagnosis)
  Infil.dat$Group <- pat.new$pat.diagnosis
  ##ggplot
  ggplot(infil.melt,aes(x=Immune.cell, y=Estimation,fill=Group))+
    geom_boxplot(outlier.size = 0.5)+
    ggtitle(title)+ 
    scale_fill_manual(values = c("#377EB8",brewer.pal(n = 9, name = "YlOrRd")[c(4,6,8,9)]))+
    #annotate("text",x = 1:n, y = as.numeric(apply(Infil.est[2:(n+1)],2,max))+0.1, label = p.star)+
    theme_bw()+
    ylab("Immune cell Abundance")+
    xlab("Immune Cell")+
    theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
          axis.text.x=element_text(angle=30,size=xsize,face = "bold",hjust=1),
          axis.text.y=element_text(size=xsize,face = "bold",hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12,face = "bold"),
          legend.text = element_text(size = 12),
          legend.position = "bottom")
  
}

##function of BCR clustering 
SeqDist <- function(x, y) {
  if (nchar(x) != nchar(y))
    return(NA)
  x.l = unlist(strsplit(x, ''))
  y.l = unlist(strsplit(y, ''))
  nn1 = length(which(x.l != y.l))
  nn2 = length(which(x.l == '-' & y.l != '-'))
  nn3 = length(which(x.l != '-' & y.l == '-'))
  return(nn1 - nn2 - nn3)
}
SeqDist.AA <- function(x, y) {
  if (nchar(x) != nchar(y))
    return(NA)
  x.l = unlist(strsplit(x, ''))
  y.l = unlist(strsplit(y, ''))
  tmp.vv = which(x.l != '-' & y.l != '-')
  tmp.st.x = min(which(x.l != '-'))
  tmp.ed.x = max(which(x.l != '-'))
  tmp.st.y = min(which(y.l != '-'))
  tmp.ed.y = max(which(y.l != '-'))
  gap.vv = which(x.l == '-' &
                   y.l != '-' | x.l != '-' & y.l == '-')## needs work
  gap.vv = gap.vv[gap.vv >= max(c(tmp.st.x, tmp.st.y)) &
                    gap.vv <= min(c(tmp.ed.x, tmp.ed.y))]
  tmp = sum(diag(BLOSUM62[x.l[tmp.vv], y.l[tmp.vv]]))
  tmp0 = sum(diag(BLOSUM62[x.l[tmp.vv], x.l[tmp.vv]]))
  score = (tmp0 - tmp + length(gap.vv)) / length(tmp.vv)
  return(score)
}
MergeSeqs <- function(seqs, dd) {
  
  tmp.seqs = gsub('-', '', seqs)
  tmp.vv = order(nchar(tmp.seqs))
  seqs = seqs[tmp.vv]
  dd = dd[tmp.vv, ]
  tmp.seqs = tmp.seqs[tmp.vv]
  newseqs = tmp.seqs
  seq.labs = rep(0, length(newseqs))
  nn = length(seq.labs)
  while (T) {
    for (ii in 1:nn) {
      if (length(grep(newseqs[ii], newseqs[-ii])) > 0)
        seq.labs[ii] = 1
    }
    if (sum(seq.labs) == 0)
      break
    vv0 = which(seq.labs == 0)
    newseqs = newseqs[vv0]
    seq.labs = seq.labs[vv0]
    nn = length(newseqs)
    seqs = seqs[vv0]
    dd = dd[vv0, ]
  }
  return(list(SS = seqs, DD = dd))
}
CreateMotifList <- function(mm) {
  tmp.mat = matrix(unlist(strsplit(rep(mm, 8, ''), '')), nrow = 8, byrow =
                     T)
  diag(tmp.mat) = '.'
  mm.list = apply(tmp.mat, 1, paste, collapse = '')
  return(mm.list)
}
MergeMotifs <- function(motif.list) {
  ## Merge motifs by allowing one mismatch
  unique.motifs = c()
  for (mm in motif.list) {
    mm.list = CreateMotifList(mm)
    sign = 0
    for (tmp.mm in mm.list) {
      if (length(grep(tmp.mm, unique.motifs)) > 0) {
        sign = 1
        break
      }
    }
    if (sign == 0)
      unique.motifs = c(unique.motifs, mm)
  }
  return(unique.motifs)
}

BuildBCRlineage <- function(sampleID, Bdata = BCRdata, start=3, end=10) {
  ## Given sample ID, start and end position of complete CDR3,  return all the lineages in the sample
  # sampleID <- "SRR3184301"
  # Bdata = cdr3.bcr.heavy
  # start=3
  # end=10
  tmp.dd.ss = subset(Bdata, sample == sampleID)
  tmp.dd.ss = tmp.dd.ss[!duplicated(tmp.dd.ss[,"CDR3nt"]),]
  if (is.null(dim(tmp.dd.ss)))
    return(NA)
  tmp.comp.vv <- which(tmp.dd.ss[, "is_complete"] == "Y")
  comp.CDR3.ss = data.frame(CDR3aa = tmp.dd.ss[tmp.comp.vv, "CDR3aa"])
  if (length(comp.CDR3.ss) == 0)
    return(NA)
  tmp.tt = table(substr(comp.CDR3.ss$CDR3aa, start, end))
  tmp.tt = sort(tmp.tt, decreasing = T)
  tmp.tt <- tmp.tt[which(nchar(names(tmp.tt))==(end-start+1))]
  tmp.motifs = MergeMotifs(names(tmp.tt))
  count = 0
  BCRlineage = c() ## a list of BCR lineage trees
  kept.motifs = c()
  for (mm in tmp.motifs) {
    mm.list = CreateMotifList(mm)
    tmp.vv.ss = c()
    for (tmp.mm in mm.list) {
      tmp.vv.ss = c(tmp.vv.ss, grep(tmp.mm, tmp.dd.ss$CDR3aa))
    }
    tmp.vv.ss = unique(tmp.vv.ss)
    if (length(tmp.vv.ss) < 2)
      next
    SEQs = unique(as.character(tmp.dd.ss[tmp.vv.ss, "CDR3nt"]))
    #SEQs = SEQs$CDR3nt
    tmp.dd0 = tmp.dd.ss[tmp.vv.ss, ]
    setDF(tmp.dd0)   ###format as dataframe
    rownames(tmp.dd0) = tmp.dd0$CDR3nt   ####same cdr3dna, same cdr3aa, different Ig gene and totaldna
    tmp.dd0 = tmp.dd0[SEQs, ]
    if (length(SEQs) < 3)
      next
    MSAalign = msa(DNAStringSet(SEQs), 'ClustalW')
    seqs = as.character(attributes(MSAalign)$unmasked)
    seqs0 = gsub('-', '', seqs)
    tmp.dd0 = tmp.dd0[match(seqs0, SEQs),]
    tmp = MergeSeqs(seqs, tmp.dd0)
    seqs = tmp$SS
    tmp.dd0 = tmp$DD
    if (is.null(dim(tmp.dd0)))
      next
    nn = nrow(tmp.dd0)
    if (nn <= 3)
      next
    sDist = matrix(0, nn, nn)
    for (ii in 1:nn) {
      for (jj in ii:nn) {
        if (jj == ii)
          next
        tmp.dist = SeqDist(seqs[ii], seqs[jj])
        sDist[ii, jj] = sDist[ii, jj] + tmp.dist
      }
    }
    kept.motifs = c(kept.motifs, mm)
    rownames(tmp.dd0) = NULL
    lineage.obj = list(distMat = sDist,
                       Sequences = seqs,
                       data = tmp.dd0)
    BCRlineage = c(BCRlineage, list(lineage.obj))
    count = count + 1
  }
  names(BCRlineage) = kept.motifs
  return(BCRlineage)
}

### BCR cluster & isotype class switch in each sample
get.bcr.cluster.classswitch <- function(bcr_clusters){
  bcr.cluster.isotypes <- NULL
  for(i in 1:length(bcr_clusters)){
    tmp <- bcr_clusters[[i]]
    if(is.null(tmp)==T){
      next
    }
    if(is.na(tmp)==T){
      next
    }
    tmp.cw <- matrix(0, nrow=length(tmp), ncol=12)
    colnames(tmp.cw) <- c('filename', 'motif', 'IGHA1','IGHA2','IGHE','IGHD','IGHM','IGHG1','IGHG2','IGHG3','IGHG4', 'Unidentified')
    tmp.cw[,'filename'] <- names(bcr_clusters)[i]
    tmp.cw[,'motif'] <- names(tmp)
    for(j in 1:length(tmp)){
      # j <- 1
      tmp_cluster <- tmp[[j]]
      tmp.is <- as.character(tmp_cluster$data$C)
      id <- which(tmp.is %in% c('IGHA1','IGHA2','IGHE','IGHD','IGHM','IGHG1','IGHG2','IGHG3','IGHG4'))
      if(length(id) > 0){
        tmp.is[-id] = 'Unidentified'
      }else{
        tmp.is = rep('Unidentified', length(tmp.is))
      }
      isotype.count <- table(tmp.is)
      tmp.cw[j, names(isotype.count)]=isotype.count
    }
    bcr.cluster.isotypes <- rbind(bcr.cluster.isotypes, tmp.cw)
  }
  return(bcr.cluster.isotypes)
}


####function of compare and boxplot
TwoLevelCompare <- function(df2,name1,name2){
  print(paste("compare",name2,"level based on two groups of",name1,sep = " "))
  df <- cbind.data.frame(Expr = df2, Group = names(df2))
  p <- signif(wilcox.test(Expr ~ Group,data = df)$p.value,4)
  print(p)
  text <- as.character(cut(p, right = FALSE, breaks = c(0, 0.001, 0.01, 0.05,1), labels = c("***","**","*"," ")))
  gplot <- ggplot(df, aes(Group, Expr,fill = Group)) +
    geom_boxplot() + 
    #geom_point(size = 0.5)+
    geom_jitter(shape=16, position=position_jitter(0.2),size = 1)+
    theme_bw() +
    scale_fill_manual(values = c("#E41A1C","#377EB8"))+
    ggtitle("") +
    xlab(paste(name1,"%", sep = " "))+
    ylab(paste(name2," Infiltration",sep = ""))+
    annotate("text",label = text,x = 1.5, y = max(df$Expr)*0.95,size = 7,color = "red")+
    theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
          axis.text.x = element_text(size = 12,face = "bold", hjust = 0.5, angle = 0),
          axis.text.y = element_text(size = 12,face = "bold"),
          axis.title.x = element_text(size = 12,face = "bold"),
          axis.title.y = element_text(size = 12,face = "bold"),
          legend.position = "none") 
  return(gplot)
}

###function of extrating different rows by comparing two dataframes
uniqdf <- function(df1,df2){
  df1[!duplicated(rbind(df2, df1))[-(1:nrow(df2))],]
}

##fisher exact test for mutated genes
MutGeneFisher <- function(row_split,g,cl1,cl2){
  cl1.pats <- rownames(row_split)[which(row_split$Cluster == cl1)]
  sub.mut1 <- as.character(unlist(oncomatrix.output[grep(g,rownames(oncomatrix.output),value = TRUE),cl1.pats]))
  g.mut1 <- length(which(sub.mut1 != ""))
  cl2.pats <- rownames(row_split)[which(row_split$Cluster == cl2)]
  sub.mut2 <- as.character(unlist(oncomatrix.output[grep(g,rownames(oncomatrix.output),value = TRUE),cl2.pats]))
  g.mut2 <- length(which(sub.mut2 != ""))
  a <- g.mut1
  b <- length(sub.mut1) - a
  c <- g.mut2
  d <- length(sub.mut2) - c
  fisher.p <- fisher.test(matrix(c(a,b,c,d),nrow = 2))$p.value
  logOddsRatio <- log2(a*d+1) - log2(b*c+1)
  res <- cbind.data.frame(gene = g, cl1 = cl1, cl2 = cl2, fisher.p = fisher.p, logOddsRatio = logOddsRatio )
  return(res)
}

##compare data between clusters
SplitComparePlot <- function(row_split,idx,feature,y_pos){
  split.alt <- row_split %>% mutate(Index = idx)
  split.alt$Cluster <- factor(split.alt$Cluster, levels = c("C1", "C2", "C3"))
  gp <- ggplot(split.alt, aes(x = Cluster, y = Index, fill = Cluster))+
    geom_violin(trim = TRUE,alpha = 0.6)+
    geom_boxplot(fill = "white")+
    #geom_point(position = position_jitterdodge())+
    theme_bw() +
    ylab("Score")+
    ggtitle(feature)+
    scale_fill_manual(values = brewer.pal(n = 3, name = "Set3"))+
    stat_compare_means(comparisons = list(c("C1", "C2"), c("C2", "C3"), c("C1", "C3")),  #, c("Cluster2", "Cluster3"), c("Cluster1", "Cluster3")
                       method = "wilcox.test",
                       label = "p.signif",
                       tip.length = 0.01,
                       # step_increase = 0.05,
                       color = "red",
                       size = 5,
                       hide.ns = TRUE,
                       #label.y = y_pos,
                       label.y.npc = "bottom",
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c( "***", "**", "*", "ns")))+
    stat_compare_means(label.y = y_pos, label.x = 1.5, size = 3)+
    theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
          axis.text.x=element_text(angle=70,size=12,face = "bold",hjust=1),
          axis.text.y=element_text(size=12,face = "bold",hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12, face = "bold"),
          legend.position = "none")
  print(gp)
}

###function of copmuting the correlation between gene expression(weighted by expression multiply rank) and immune infiltration corrected by tumor purity
CorrelatinExprInfilCorrectPurity <- function(cells,g.list,Expr,Infiltration,purity,method){
  # ge.select <- apply(Expr,2,function(x)match(x,sort(x))) #/dim(ge.rank.zscore)[2]  ###rank of genes
  # rownames(ge.select) <- rownames(Expr)
  CorrCellGene <- lapply(cells, function(c){
    CorrGes <- lapply(g.list, function(g){
      if(length(g) > 1 & method == "/"){
        ge.mul <- Expr[g,]#*ge.select[g,]
        ge.com <- log2((ge.mul[2,]+1)/(ge.mul[1,]+1))
        #ge.mul <- apply(ge.select[g,], 1, function(x) log2(x[[2]]*Expr[g[2],]/x[[1]]*Expr[g[1],]+1))  
        name <- paste(rev(g),collapse = "/")
      }
      if(length(g) > 1 & method == "+"){
        #ge.mul <- ge.select[g,]*Expr[g,]
        ge.mul <- Expr[g,] #*ge.select[g,]
        ge.com <- log2(colMeans(ge.mul)+1)
        #ge.mul <- apply(ge.select[g,], 1, function(x)log2(x[[2]]*Expr[g[2],]+x[[1]]*Expr[g[1],]+1))  ##add up rank
        name <- paste(g,collapse = "+")
      }
      if(length(g) == 1 & method == "+"){
        ge.mul <- Expr[g,] #*ge.select[g,]
        ge.com <- unlist(log2(ge.mul +1))
        name <- g
      }
      # if(length(g) == 1 & method == ""){
      #   ge.mul <- ge.select[g,]*Expr[g,]
      #   ge.com <- unlist(log2(ge.mul +1))
      #   name <- g
      # }
      names(ge.com) <- colnames(Expr)
      #ge.rank <- sort(ge.select)
      dat <- cbind.data.frame(Expression = (unlist(ge.com) - min(unlist(ge.com)))/(max(unlist(ge.com)) - min(unlist(ge.com))), #unlist(ge.com),
                              # Expression.zscore =  (unlist(ge.com) - mean(unlist(ge.com)))/sd(unlist(ge.com)),
                              Infil = unlist(Infiltration[c,names(ge.com)]),
                              #Infil.zscore = (unlist(Infiltration[c,names(ge.com)]) - mean(unlist(Infiltration[c,names(ge.com)])))/sd(unlist(Infiltration[c,names(ge.com)])),
                              purity = purity[names(ge.com)],
                              Cell = c)
      test <- pcor.test(dat$Expression,dat$Infil,dat$purity, method = "spearman")
      Correlation = signif(test$estimate,2)
      Pval <- test$p.value
      p.adj <- p.adjust(Pval, method = "BH")
      pstar = as.character(cut(p.adj, right = FALSE, breaks = c(0, 0.001, 0.01, 0.05,1), labels = c("***","**","*"," ")))
      res <- cbind.data.frame(Gene = name, correlation = Correlation, padj = p.adj, p = pstar, cell = c)
      return (res)
    })
    EachCell <- do.call("rbind",CorrGes)
    return(EachCell)
  })
  all.corr <- do.call("rbind",CorrCellGene)
  all.corr$p <- as.character(all.corr$p)
  return(all.corr)
}


###function of comparing TCR and BCR metrics
CompareTBCRmetric <- function(dat, feature){
  ###correlation
  cor.test <- cor.test(dat$F1,dat$F2,method = "spearman")
  cor <- signif(cor.test$estimate,2)
  p <- cut(cor.test$p.value,breaks = c(0, 0.001, 0.01, 0.05,1), labels = c("***","**","*"," "))
  ###scatter plot
  ggplot(dat,aes(x = F1, y = F2))+
    geom_point()+
    geom_smooth(method = "lm")+
    theme_bw()+
    xlab("TCR")+
    ylab("BCR")+
    # scale_color_manual(values = c("#E41A1C", "#377EB8"), name = "")+
    # coord_cartesian(xlim = 1.3 * c(min(dat$F2), max(dat$F2)), 
    #                 ylim = 1.3 * c(min(dat$F1), max(dat$F1))) +   
    # geom_encircle(data = dat %>% dplyr::filter(F3 == "Immune_High"), aes(x = F2, y = F1))+
    # geom_encircle(data = dat %>% dplyr::filter(F3 == "Immune_Low"), aes(x = F2, y = F1))+
    ggtitle(paste(feature,paste(paste(paste("Correlation = ", cor, sep = ""), p, sep = "  ")),sep = "\n"))+
    theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
          axis.text.x=element_text(angle=70,size=12,face = "bold",hjust=1),
          axis.text.y=element_text(size=12,face = "bold",hjust=1),
          axis.title.x = element_text(size=12,face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold")
    ) 
}

###function of comparing metrics in patients from low to high immune score
DotLowHigh <- function(pat.ord, feature, score.aver, idx){
  ##given a rank of patients
  ##extract metric
  df <- cbind.data.frame(pat = pat.ord, idex = idx, 
                         group = names(pat.ord),isc = score.aver[pat.ord])
  df$group <- as.character(df$group)
  df.cl.ord<- df %>% arrange(isc)
  # ##given a rank of groups
  # gr.ord <- df.ord %>% group_by(group) %>% dplyr::summarise(gr.levl = mean(idex)) %>% arrange(gr.levl)
  # ##rank patients in each cluster
  # df.cl.ord <- NULL
  # for(c in gr.ord$group){
  #   tmp <- subset(df.ord, group == c)
  #   df.cl.ord <- rbind(df.cl.ord,tmp)
  # }
  # ###scatter plot
  ###rect data
  md.idex1 <- min(subset(df.cl.ord, isc >= median(df.cl.ord$isc), select = "idex"))
  md.idex2 <- max(subset(df.cl.ord, isc <= median(df.cl.ord$isc), select = "idex"))
  df.rect <- cbind.data.frame(x1 = c(min(df.cl.ord$isc),median(df.cl.ord$isc)),
                              x2 = c(median(df.cl.ord$isc),max(df.cl.ord$isc)),
                              y1 = c(min(df.cl.ord$idex),md.idex1),
                              y2 = c(md.idex2,max(df.cl.ord$idex)),
                              group = c("Immune_Low","Immune_High"))
  ###calculate correlation
  ct <- cor.test(df.cl.ord$idex,df.cl.ord$isc, method = "spearman")
  ct.cor <- signif(ct$estimate,2)
  ct.p <- as.character(cut(signif(ct$p.value,2),breaks = c(0, 0.001, 0.01, 0.05,1), labels = c("***","**","*"," "),include.lowest = TRUE))
  gp <- ggplot(df.cl.ord) +
    # geom_line()+
    #geom_rect(data = df.rect, mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill= group),alpha = 0.3)+
    geom_vline(xintercept = median(df.cl.ord$isc), linetype = "dashed", color = "black")+
    geom_point(aes(x=isc, y=idex, color = group), alpha = 0.8)+
    geom_smooth(aes(x=isc, y=idex), method = lm, color = "black")+
    scale_color_manual(values = c("#E41A1C","#377EB8"))+
    ggtitle(paste(feature, paste(paste("Corr = ", ct.cor, sep = ""), ct.p, sep = " "), sep = "\n"))+
    xlab(expression(Immune_Low %->% Immune_High))+
    ylab("Score")+
    theme_bw()+
    # annotate("text",x = (min(df.cl.ord$isc-median(df.cl.ord$isc)))/2, y = (md.idex2+max(df.cl.ord$idex))/2,
    #          label = paste(paste("Corr = ", ct.cor, sep = ""), ct.p, sep = "\n"), size =5) +
    theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=10,r=10,t=10,b=10),face = "bold", colour = "black"),
          axis.text.x=element_text(color = "white"),
          axis.text.y=element_text(size=12,face = "bold",hjust=1),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.ticks.x = element_blank(),
          legend.position = "None")
  mp <- ggboxplot(df.cl.ord, x = "group", y = "idex",color = "group", palette = c("#377EB8","#E41A1C"),add = "jitter")+
    stat_compare_means(comparisons = list(c("Immune_High", "Immune_Low")), method = "wilcox.test", label.y = max(df.cl.ord$idex)-0.05,
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c( "***", "**", "*", "ns")))+
    xlab("")+
    ggtitle("")+
    #theme_light()+
    theme(legend.position = "None",
          plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=10,r=10,t=10,b=10),face = "bold", colour = "black"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(colour='white'),
          axis.text.x = element_text(colour='white'),
          axis.ticks.x=element_blank(),
          axis.ticks.y = element_blank())
  cp <- ggarrange(gp, mp,ncol=2, nrow=1, widths = c(4,1),heights = c(4,4))
  return(cp)
}

###function of comparing metrics in patients from low to high immune score
DotLowHighRemoveBox <- function(pat.ord, feature, score.aver, idx){
  ##given a rank of patients
  ##extract metric
  df <- cbind.data.frame(pat = pat.ord, idex = idx, 
                         group = names(pat.ord),isc = score.aver[pat.ord])
  df$group <- as.character(df$group)
  df.cl.ord<- df %>% arrange(isc)
  ###calculate correlation
  ct <- cor.test(df.cl.ord$idex,df.cl.ord$isc, method = "spearman")
  ct.cor <- signif(ct$estimate,2)
  ct.p <- as.character(cut(signif(ct$p.value,2),breaks = c(0, 0.001, 0.01, 0.05,1), labels = c("***","**","*"," "),include.lowest = TRUE))
  gp <- ggplot(df.cl.ord) +
    # geom_line()+
    #geom_rect(data = df.rect, mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill= group),alpha = 0.3)+
    geom_vline(xintercept = median(df.cl.ord$isc), linetype = "dashed", color = "black")+
    geom_point(aes(x=isc, y=idex, color = group), alpha = 0.8)+
    geom_smooth(aes(x=isc, y=idex), method = lm, color = "black")+
    scale_color_manual(values = c("#E41A1C","#377EB8"))+
    ggtitle(paste(feature, paste(paste("Corr = ", ct.cor, sep = ""), ct.p, sep = " "), sep = "\n"))+
    xlab(expression(Immune_Low %->% Immune_High))+
    ylab("Score")+
    theme_bw()+
    # annotate("text",x = (min(df.cl.ord$isc-median(df.cl.ord$isc)))/2, y = (md.idex2+max(df.cl.ord$idex))/2,
    #          label = paste(paste("Corr = ", ct.cor, sep = ""), ct.p, sep = "\n"), size =5) +
    theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=10,r=10,t=10,b=10),face = "bold", colour = "black"),
          axis.text.x=element_text(color = "white"),
          axis.text.y=element_text(size=12,face = "bold",hjust=1),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.ticks.x = element_blank(),
          legend.position = "None")
  return(gp)
}

###function of generating oncomatrix from maf file
createOncoMatrix = function(m, g = NULL, chatty = TRUE, add_missing = FALSE){
  
  if(is.null(g)){
    stop("Please provde atleast two genes!")
  }
  
  subMaf = subsetMaf(maf = m, genes = g, includeSyn = FALSE, mafObj = FALSE)
  
  if(nrow(subMaf) == 0){
    if(add_missing){
      numericMatrix = matrix(data = 0, nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(numericMatrix) = g
      colnames(numericMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])
      
      oncoMatrix = matrix(data = "", nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(oncoMatrix) = g
      colnames(oncoMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])
      
      vc = c("")
      names(vc) = 0
      
      return(list(oncoMatrix = oncoMatrix, numericMatrix = numericMatrix, vc = vc))
    }else{
      return(NULL)
    }
  }
  
  if(add_missing){
    subMaf[, Hugo_Symbol := factor(x = Hugo_Symbol, levels = g)]
  }
  
  oncomat = data.table::dcast(data = subMaf[,.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
                              fun.aggregate = function(x){
                                x = unique(as.character(x))
                                xad = x[x %in% c('Amp', 'Del')]
                                xvc = x[!x %in% c('Amp', 'Del')]
                                
                                if(length(xvc)>0){
                                  xvc = ifelse(test = length(xvc) > 1, yes = 'Multi_Hit', no = xvc)
                                }
                                
                                x = ifelse(test = length(xad) > 0, yes = paste(xad, xvc, sep = ';'), no = xvc)
                                x = gsub(pattern = ';$', replacement = '', x = x)
                                x = gsub(pattern = '^;', replacement = '', x = x)
                                return(x)
                              } , value.var = 'Variant_Classification', fill = '', drop = FALSE)
  
  #convert to matrix
  data.table::setDF(oncomat)
  rownames(oncomat) = oncomat$Hugo_Symbol
  oncomat = as.matrix(oncomat[,-1, drop = FALSE])
  
  variant.classes = as.character(unique(subMaf[,Variant_Classification]))
  variant.classes = c('',variant.classes, 'Multi_Hit')
  names(variant.classes) = 0:(length(variant.classes)-1)
  
  #Complex variant classes will be assigned a single integer.
  vc.onc = unique(unlist(apply(oncomat, 2, unique)))
  vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
  names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)])+1), length(vc.onc))
  variant.classes2 = c(variant.classes, vc.onc)
  
  oncomat.copy <- oncomat
  #Make a numeric coded matrix
  for(i in 1:length(variant.classes2)){
    oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
  }
  
  #If maf has only one gene
  if(nrow(oncomat) == 1){
    mdf  = t(matrix(as.numeric(oncomat)))
    rownames(mdf) = rownames(oncomat)
    colnames(mdf) = colnames(oncomat)
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
  
  #convert from character to numeric
  mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
  rownames(mdf) = rownames(oncomat.copy)
  
  
  #If MAF file contains a single sample, simple sorting is enuf.
  if(ncol(mdf) == 1){
    sampleId = colnames(mdf)
    mdf = as.matrix(mdf[order(mdf, decreasing = TRUE),])
    colnames(mdf) = sampleId
    
    oncomat.copy = as.matrix(oncomat.copy[rownames(mdf),])
    colnames(oncomat.copy) = sampleId
    
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  } else{
    #Sort by rows as well columns if >1 samples present in MAF
    #Add total variants per gene
    mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
      length(x[x != "0"])
    }))
    #Sort by total variants
    mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
    #colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
    nMut = mdf[, ncol(mdf)]
    
    mdf = mdf[, -ncol(mdf)]
    
    mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix
    
    mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
    tmdf = t(mdf) #transposematrix
    mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort
    
    mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
    mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
    mdf = mdf.temp.copy
    
    #organise original character matrix into sorted matrix
    oncomat.copy <- oncomat.copy[,colnames(mdf)]
    oncomat.copy <- oncomat.copy[rownames(mdf),]
    
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
}