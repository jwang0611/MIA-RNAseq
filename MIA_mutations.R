###read in maf file
n <- 59
lung.maf <- read.maf(maf="~/Documents/Description_study/vcf2maf_data/lung_merged_remove_noveldbSNP_Silent_Flank_IGR_Intron_RNA_UTR_removeTBCR.varscan.base.snp.indel.filter.maf",
                     useAll = TRUE, isTCGA = FALSE,removeDuplicatedVariants = TRUE)
col = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')

####mutation overview of MIA
pdf(file = "~/Documents/Description_study/figures_v2/figures2_mutation_summary_plot.pdf",height = 4,width = 6)
plotmafSummary(maf = lung.maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE, color = col,textSize = 3,titleSize = c(1,1))
dev.off()

####top mutated genes of MIA
pdf(file = "~/Documents/Description_study/figures_v2/figures2_MIA_oncoplot2.pdf", height = 10, width = 9)
oncoplot(maf = lung.maf, top = 30, fontSize = 16,colors = col, removeNonMutated = FALSE,
         showTumorSampleBarcodes = FALSE, drawColBar = FALSE) #sampleOrder = pat.order, 
dev.off()

####transition and tranversion rate
lung.titv = titv(maf = lung.maf, plot = FALSE, useSyn = TRUE)
#plot titv summary
pdf(file = "~/Documents/Description_study/figures_v2/figures2_MIA_titv.pdf", height = 4, width = 4)
plotTiTv(res = lung.titv)
dev.off()

####mutation of EGFR, KRAS and TP53
pdf(file = "~/Documents/Description_study/figures_v2/figures2_MIA_KRAS_EGFR_TP53.pdf", height = 3, width = 6)
oncoplot(maf = lung.maf,genes= c("TP53","KRAS","EGFR"), removeNonMutated = FALSE,colors = col,showTumorSampleBarcodes = FALSE)
dev.off()


##############------check msi score and CT transition -----------------------------#####
######################################################################################################################################
###read in msi score
msi <- read.table("~/Documents/Description_study/data/MIA_Pat_msi_score.txt", sep = "\t")
###titv
titv = titv(maf = lung.maf, useSyn = TRUE, plot = FALSE)
titv.counts = data.frame(titv$raw.counts)%>% 
    mutate(total.titv = rowSums(.[2:7])) %>%
    mutate(Patient = paste("Pat",as.character(patient.info[match(Tumor_Sample_Barcode, patient.info$tumor),"patientID"]), sep = "")) 
colnames(titv.counts)[2:7] <- gsub(">","",colnames(titv.counts)[2:7])
titv.counts <- titv.counts %>% mutate(CT.rate = (C.T + T.C)/total.titv) %>% 
    mutate(msi.score = msi[match(Patient, msi$V1),"V2"]) %>%
    mutate(group = ifelse(CT.rate > median(CT.rate),"High CT transition", "Low CT transition"))
head(titv.counts)
##boxplot
pdf("~/Documents/Description_study/figures_v2/msi_score_CT_transition.pdf", width = 4, height = 4)
ggplot(titv.counts, aes(x = group, y = msi.score, fill = group))+
    geom_boxplot()+
    stat_compare_means()+
    theme_bw()+
    theme(legend.position = "top",
        axis.text.x=element_text(size = 12,face = "bold", angle = 30, hjust = 1),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.ticks.x.bottom = element_blank())
dev.off()

