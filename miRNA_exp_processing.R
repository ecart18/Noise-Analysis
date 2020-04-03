
rm(list = ls())
work_directory <- getwd()


##### miRNA-seq Data processing
mESCs_miRNA_file_path <- paste(work_directory,'/data/GSE89690_miRNAseq_expression.TSV',sep = '')
mESCs_miRNA <- read.table( mESCs_miRNA_file_path,stringsAsFactors = FALSE)
COLNAMES <- c('miRNAname','R2i_RB18','R2i_RB20','Serum_RB18','Serum_RB20','2i_RB18_V2','2i_RB20_V2')
COLNAMES <- c('miRNAname','Ddx5rep1','Ddx5rep2','Ddx5line1','Ddx5line2','Ddx5line3','Ddx5line4')
colnames(mESCs_miRNA) <- COLNAMES
mESCs_miRNA <- mESCs_miRNA[-1,]
mESCs_miRNA[,2:7] <- lapply(mESCs_miRNA[,2:7],as.numeric)
meanvalue <- c((mESCs_miRNA$Serum_RB18 + mESCs_miRNA$Serum_RB20)/2)
meanvalue <- c((mESCs_miRNA$Ddx5rep1 + mESCs_miRNA$Ddx5rep2)/2)
mESCs_miRNA <- cbind(mESCs_miRNA,meanvalue)
mESCs_miRNA <- mESCs_miRNA[order(mESCs_miRNA$meanvalue,decreasing=T),]

### find miRNA-seq belongs to mmu and generate modifyname------Dataframe : mESCs_miRNA_mmu
mmuindex <- grep("mmu-*",mESCs_miRNA$miRNAname)
mESCs_miRNA_mmu <- mESCs_miRNA[mmuindex,]
mESCs_miRNA_mmu$modifyname <- mESCs_miRNA_mmu$miRNAname
modifyname <- c()
for (i in 1:length(mESCs_miRNA_mmu$miRNAname) ){
  text <- mESCs_miRNA_mmu$miRNAname[i]
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
 
mmuindex2 <- grep("[\\*]$",modifyname)
for (i in mmuindex2){
  text <- modifyname[i]
  modifyname[i] <- substr(text,1,nchar(text)-1)
}
mESCs_miRNA_mmu$modifyname <- modifyname

### unique and filt------Dataframe : mESCs_miRNA_mmu_unique & mESCs_miRNA_mmu_fiter
mESCs_miRNA_mmu <- mESCs_miRNA_mmu[order(mESCs_miRNA_mmu$meanvalue,decreasing=T),]
index<-duplicated(mESCs_miRNA_mmu$modifyname) 
mESCs_miRNA_mmu_unique <- mESCs_miRNA_mmu[!index,]
mESCs_miRNA_mmu_fiter <- mESCs_miRNA_mmu_unique[mESCs_miRNA_mmu_unique$meanvalue >= 0,]
save(mESCs_miRNA,mESCs_miRNA_mmu,mESCs_miRNA_mmu_unique, mESCs_miRNA_mmu_fiter, file=paste(work_directory,'/data/mESCs_miRNA_GSE89690.RData',sep = ''))

