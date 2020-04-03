rm(list = ls())
work_directory <- getwd()

source(paste(work_directory,'/utils/noise_decomposition_function.R',sep = ''))
source(paste(work_directory,'/utils/utils.R',sep = ''))

detachAllPackages()

param_mESCs_scRNA_filename <-  paste(work_directory, '/data/mESCs_scRNA_islam_deconv.Rdata', sep = '')
param_miRNAfilename <-  paste(work_directory, '/data/mESCs_miRNA_GSE89690.Rdata', sep = '')
param_savefile <-  paste(work_directory, '/data/mESCs_miRNA_TargetScan_mapping.Rdata ', sep = '')
param_sitetype <- 'Allmer'
param_rankpercent <- 101


load(param_mESCs_scRNA_filename)
load(param_miRNAfilename)
mESCs_miRScan_file_path <- paste(work_directory,'/data/Conserved_Site_Context_Scores.txt',sep = '')
mESCs_miRScan <- read.table( mESCs_miRScan_file_path ,skip = 1, stringsAsFactors = FALSE)
colnames(mESCs_miRScan) <- c('GeneID','GeneSymbol','TranscriptID','GeneTaxID',
                             'miRNA','SiteType','UTRstart','UTRend','contextscore',
                             'contextscorepercentile','weightedcontextscore','weightedcontextscorepercentile')
mESCs_miRScan[,6:12] <- lapply(mESCs_miRScan[,6:12], as.numeric)
mESCs_miRScan <- mESCs_miRScan[!is.na(mESCs_miRScan$weightedcontextscorepercentile),]
mESCs_miRScan_file_path2 <- paste(work_directory,'/data/Nonconserved_Site_Context_Scores.txt',sep = '')
mESCs_miRScan2 <- read.table(mESCs_miRScan_file_path2 ,skip = 1, stringsAsFactors = FALSE)
colnames(mESCs_miRScan2) <- c('GeneID','GeneSymbol','TranscriptID','GeneTaxID',
                             'miRNA','SiteType','UTRstart','UTRend','contextscore',
                             'contextscorepercentile','weightedcontextscore','weightedcontextscorepercentile')

mESCs_miRScan2[,6:12] <- lapply(mESCs_miRScan2[,6:12], as.numeric)
mESCs_miRScan2 <- mESCs_miRScan2[!is.na(mESCs_miRScan2$weightedcontextscorepercentile),]
mESCs_miRScan2 <- mESCs_miRScan2[mESCs_miRScan2$weightedcontextscorepercentile>=param_rankpercent,]

mESCs_miRScan <- rbind(mESCs_miRScan,mESCs_miRScan2)
merNum_threshod <- 6
sitetypenum <- 3  #3 for 8 mer ;   1/2 for 7 mer;   4 for 6 mer
mers <- abs(mESCs_miRScan$UTRstart - mESCs_miRScan$UTRend) + 1
mESCs_miRScan <- cbind(mESCs_miRScan,mers)

if (param_sitetype == '6mer'){
	mESCs_miRScan <- mESCs_miRScan[mESCs_miRScan$SiteType==4,]}
if (param_sitetype == '7mer'){
	mESCs_miRScan <- mESCs_miRScan[mESCs_miRScan$SiteType<sitetypenum,]}
if (param_sitetype == '8mer'){
	mESCs_miRScan <- mESCs_miRScan[mESCs_miRScan$SiteType==sitetypenum,]}
if (param_sitetype == 'Allmer'){
	mESCs_miRScan <- mESCs_miRScan}

mmu_index <- grep("mmu-*",mESCs_miRScan$miRNA)
mESCs_miRScan_mmu <- mESCs_miRScan[mmu_index,]
mESCs_miRScan_mmu$modifyname <- mESCs_miRScan_mmu$miRNA
modifyname <- c()

for (i in 1:length(mESCs_miRScan_mmu$miRNA) ){
  text <- mESCs_miRScan_mmu$miRNA[i]
  text <-  as.character(text)
  indextmp <- gregexpr("-",text)
  indextmp2 <- indextmp[[1]][3]
  if (is.na(indextmp2)){
    modifyname[i] <- as.character(text)
  }
  else{
    modifyname[i] <- substr(text,1,indextmp2-1)
  }
}

mmu_index2 <- grep("[\\*]$",modifyname)
for (i in mmu_index2){
  text <- modifyname[i]
  modifyname[i] <- substr(text,1,nchar(text)-1)
}
mESCs_miRScan_mmu$modifyname <- modifyname
GeneIDtmp <- mESCs_miRScan_mmu$GeneID
mmu_index3 <- grep("[\\.][0-9]$",GeneIDtmp)
for (i in mmu_index3){
  text <- GeneIDtmp[i]
  GeneIDtmp[i] <- substr(text,1,nchar(text)-2)
}
mESCs_miRScan_mmu$GeneID <- GeneIDtmp
dim(mESCs_miRScan_mmu)

GeneIDtmp <- mESCs_miRScan_mmu$GeneID
mmu_index3 <- grep("[\\.][0-9][0-9]$",GeneIDtmp)
for (i in mmu_index3){
  text <- GeneIDtmp[i]
  GeneIDtmp[i] <- substr(text,1,nchar(text)-3)
}
mESCs_miRScan_mmu$GeneID <- GeneIDtmp


##miRNA * mRNA matrix
df_miRNA_scRNA <- data.frame(matrix(0,nrow = nrow(mESCs_miRNA_mmu_fiter),ncol = nrow(mESCs_scRNA_removegene)))
colnames(df_miRNA_scRNA) <- rownames(mESCs_scRNA_removegene)
rownames(df_miRNA_scRNA) <- mESCs_miRNA_mmu_fiter$modifyname
df_miRNA_scRNA_value <- df_miRNA_scRNA
df_miRNA_scRNA_siteNum <- df_miRNA_scRNA
dim(df_miRNA_scRNA)
dim(df_miRNA_scRNA_value)

# overlap1  the target gene number
overlap1 <- 0
mESCs_scRNA_non_target <- mESCs_scRNA_removegene
for (i in 1:length(mESCs_miRNA_mmu_fiter$modifyname) ){
  miR_i <- mESCs_miRNA_mmu_fiter$modifyname[i]
  miR_i_index <- grep(paste('^',miR_i,'$',sep = ''),mESCs_miRScan_mmu$modifyname)
  if (length(miR_i_index) == 0) {next()}
  df_miR_i <- mESCs_miRScan_mmu[miR_i_index,]
  index<-duplicated(df_miR_i$GeneID) 
  df_miR_i <- df_miR_i[!index,]
  target_miR_i <-  df_miR_i$GeneID
  for (j in 1:length(target_miR_i) ){
    target_j <- target_miR_i[j]
    
    if ( target_j %in% rownames(mESCs_scRNA_non_target) ){
      loca <- which( rownames(mESCs_scRNA_non_target) == target_j ) 
      mESCs_scRNA_non_target <- mESCs_scRNA_non_target[-loca,]
      overlap1 <- overlap1 +1
    }
    
    if ( target_j %in% colnames(df_miRNA_scRNA) ){
      df_miRNA_scRNA[miR_i,as.character(target_j)] <- 1
      df_miRNA_scRNA_value[miR_i,as.character(target_j)] <- mESCs_miRNA_mmu_fiter[which(mESCs_miRNA_mmu_fiter$modifyname == miR_i),'meanvalue']
      df_miRNA_scRNA_siteNum[miR_i,as.character(target_j)] <- df_miR_i$mers[j]
      }
    
  }
}

ALL_target_name <- colnames(df_miRNA_scRNA)[colSums(df_miRNA_scRNA) > 0]
mESCs_scRNA_target <- mESCs_scRNA_removegene[ALL_target_name,]
# list for test every miRNA target gene 
miRNA_tests_all<-list()
miRNA_tests_all <- as.list(as.data.frame(t(df_miRNA_scRNA)))


save(miRNA_tests_all,mESCs_scRNA_non_target, 
     df_miRNA_scRNA, mESCs_scRNA_target,
     df_miRNA_scRNA_value,df_miRNA_scRNA_siteNum,
     mESCs_miRScan_mmu, mESCs_miRScan, file=param_savefile)








