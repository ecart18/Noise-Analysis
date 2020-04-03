rm(list = ls())
work_directory <- getwd()

source(paste(work_directory,'/utils/noise_decomposition_function.R',sep = ''))
source(paste(work_directory,'/utils/utils.R',sep = ''))

detachAllPackages()
library(plyr)
library(minpack.lm)
library(coin)
library(ggplot2)
library(RColorBrewer)

# noise bias test
load(paste(work_directory,'/data/mESCs_scRNA_islam_deconv.Rdata', sep=''))
load(paste(work_directory,'/data/mESCs_miRNA_GSE89690.Rdata', sep=''))
load(paste(work_directory,'/data/mESCs_miRNA_TargetScan_mapping.Rdata', sep=''))

backround = 'NonTargetCtrl'
length(miRNA_tests_all[[1]])
Result_tests_all_islam<-list()
# test for NonTargetCtrl
if(backround == 'NonTargetCtrl'){
  Result_tests_all_islam[['Mouse']]<-genesets_biased_NonTargetCtrl(miRNA_tests_all, gene_cellcycle_rank[!is.na(gene_cellcycle_rank)], df_miRNA_scRNA,mESCs_scRNA_non_target)
}
Result_tests_all_islam<-do.call('rbind', Result_tests_all_islam)
rownames(Result_tests_all_islam) <- names(miRNA_tests_all)
Result_tests_all_islam$miRexp <- mESCs_miRNA_mmu_fiter[!is.na(match(mESCs_miRNA_mmu_fiter$modifyname,names(miRNA_tests_all))),'meanvalue']
dim(Result_tests_all_islam)
Result_tests_all_islam <- Result_tests_all_islam[!is.na(Result_tests_all_islam$effect_size),]
dim(Result_tests_all_islam)
sum(rowSums(df_miRNA_scRNA)>0)
Result_tests_effects_low_islam <- Result_tests_all_islam[Result_tests_all_islam$effect_size <= 0.5,]
Islam_100_all_times1_df_miRNA_scRNA <- df_miRNA_scRNA
Islam_gene_cellcycle_rank <- gene_cellcycle_rank
Islam_mESCs_scRNA_non_target <- mESCs_scRNA_non_target
Islam_mESCs_miRNA_mmu_fiter <- mESCs_miRNA_mmu_fiter
Islam_mESCs_scRNA_removegene <- mESCs_scRNA_removegene



# adjusted effect size
negnum <- dim(mESCs_scRNA_non_target)[1]
posnum <- Result_tests_all_islam$NumTarget
sd <- sqrt( (((negnum+1)/posnum) + 1)/12/negnum )
Result_tests_all_islam$adj_effect_size <- (0.5-Result_tests_all_islam$effect_size)/sd

Islam_100_all_pvalue <- rownames(Result_tests_effects_low_islam[Result_tests_effects_low_islam$p_value<0.05,])
length(Islam_100_all_pvalue)
Islam_100_all_fdr <- rownames(Result_tests_effects_low_islam[Result_tests_effects_low_islam$FDR <= 0.1,])
length(Islam_100_all_fdr)

ttmp <- as.data.frame(matrix(mESCs_scRNA_removegene$predictedCount, nrow=dim(Islam_100_all_times1_df_miRNA_scRNA)[1], ncol=length(mESCs_scRNA_removegene$predictedCount), byrow=T))*Islam_100_all_times1_df_miRNA_scRNA
rownames(ttmp) <- rownames(Islam_100_all_times1_df_miRNA_scRNA)
Result_tests_all_islam$Sumexp <- (rowSums(ttmp)[rownames(Result_tests_all_islam)])
Result_tests_all_islam$`Dataset:` <- 'Islam et al.'
Result_tests_all_islam$miRname <- rownames(Result_tests_all_islam)
rownames(Result_tests_all_islam) <- NULL
colnames(Result_tests_all_islam)



# plot 
common_denoise_miRNA <- c("mmu-miR-455", "mmu-miR-7a", "mmu-miR-184", "mmu-miR-322", "mmu-miR-152",
                          "mmu-miR-324", "mmu-miR-15a", "mmu-miR-301b", "mmu-miR-24", "mmu-miR-15b",
                          "mmu-miR-130b", "mmu-miR-148b", "mmu-miR-301a", "mmu-miR-148b", "mmu-miR-16")
#Bubble colour
plot_df <- Result_tests_all_islam
plot_df$FDROverlap<-'FDRSpecific'
plot_df$FDROverlap[!is.na(match(plot_df$miRname,common_denoise_miRNA))]<-'FDROverlap'
plot_df$FDROverlap<-as.factor(plot_df$FDROverlap)
plot_df <- plot_df[plot_df$FDROverlap=='FDROverlap',]
#Plot index
temp<-plot_df[!duplicated(plot_df$miRname),]
temp_order<-order(temp$miRexp, temp$miRname)
plot_df$index<-factor(plot_df$miRname, levels=temp$miRname[temp_order])

d <- ggplot(plot_df, aes(index, adj_effect_size,group = interaction(index))) + 
  geom_point(color='orange', size=4,position = position_dodge(width = 0.75)) +
  geom_hline(yintercept = 0.0, linetype=2, size=0.5) + 
  scale_y_continuous(breaks = c(-2,0,2,4),limits = c(-2,5))+
  xlab('miRNA name') +
  ylab('Adjusted effect size') + 
  labs(fill="Datasets:") +
  theme_bw() + 
  theme(legend.position=c(0.85,0.10),legend.text = element_text(size = 15, hjust = 3, vjust = 3))+  
  theme(axis.title.y=element_text(angle=90, colour="black",size=15, vjust=0.5),axis.title.x=element_text(angle=0, colour="black",size=15, hjust=0.5)) +
  theme(axis.text.x = element_text(angle=90, colour="black",hjust=1,size=12))+   
  theme(axis.text.y = element_text(angle=90, colour="black",vjust=0.5,size=12)) 
d + geom_blank()
