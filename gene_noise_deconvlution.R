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

# load the estimated number of spiked-in mRNA molecules of the 92 ERCC spike-ins from ercc_counts.txt
ERCCpath <- paste(work_directory,'/data/GSE46980_ERCC_count.txt',sep = '')
erccNumber <- read.table(ERCCpath, stringsAsFactors = FALSE, skip = 1)
erccNumber <- erccNumber[,c(1,6,10)]
colnames(erccNumber) <- c('RE.sort.ID','mappingName','molecularPerChamble')
erccNumber$molecularPerChamble <- 1*erccNumber$molecularPerChamble
rownames(erccNumber) <- erccNumber[,1]
erccSeqID <- sort(rownames(erccNumber))

#load the number of sequenced mRNA molecules of the spike-ins for each cell from ercc_cell_counts.txt
erccCountpath <- paste(work_directory,'/data/GSE46980_ERCC_cell_count',sep = '')
erccCount = read.table(erccCountpath, stringsAsFactors = FALSE, skip = 1)
rownames(erccCount) = erccCount[,1]
filter_name <- match(rownames(erccCount),erccNumber$mappingName)
erccCount <- erccCount[!is.na(filter_name),]
filter_name <- filter_name[!is.na(filter_name)]
erccCount$V1 <- erccNumber$RE.sort.ID[filter_name]
rownames(erccCount) <- erccCount$V1 
erccCount <- erccCount[,8:length(erccCount[1,])]
COL_NAME <- c('A01','B01','C01','D01','E01','F01','G01','H01','A02','B02','C02','D02','E02','F02','G02','H02',	
              'A03','B03','C03','D03','E03','F03','G03','H03','A04','B04','C04','D04','E04','F04','G04','H04','A05','B05',
              'C05','D05','E05','F05','G05','H05','A06','B06','C06','D06','E06','F06','G06','H06','A07','B07','C07','D07',
              'E07','F07','G07','H07','A08','B08','C08','D08','E08','F08','G08','H08','A09','B09','C09','D09','E09','F09',
              'G09','H09','A10','B10','C10','D10','E10','F10','G10','H10','A11','B11','C11','D11','E11','F11','G11','H11',	
              'A12','B12','C12','D12','E12','F12','G12','H12')
colnames(erccCount) <- COL_NAME
erccCount[,1:length(erccCount[1,])] <- lapply(erccCount[,1:length(erccCount[1,])],as.numeric)

# quality control
index <- c("A03","B05","B06","B07","B09","B10","B12","C01","C02","C06","C08","A04",
           "C09","C10","C11","C12","D01","D05","D06","D07","D11","E07","A05","E09",
           "E12","F05","F10","G04","G05","G09","H04","H05","H06","A07","H08","H10",
           "A08","A10","A11","B01","B04")
erccCount <- erccCount[,index]
erccCount <- erccCount[rownames(erccNumber),]
erccNumber <- erccNumber[!is.na(erccCount[,1]),]
erccNumber <- erccNumber[,-2]
erccCount <- erccCount[!is.na(erccCount[,1]),]
erccSeqID <- sort(rownames(erccNumber))
colnames(erccCount) <- paste('SC_2i_',seq(1,41),sep = '')

dim(erccCount)
dim(erccNumber)
head(erccCount)
head(erccNumber)



# load genes counts
#load the number of sequenced mRNA molecules of genes for each cell from gene_cell_counts.txt.
load(paste(work_directory,'/data/allgenes_info.Rdata', sep=''))
countGenespath <- paste(work_directory, '/data/GSE46980_CombinedMoleculeCounts',sep = '')
countGenes <- read.table(countGenespath, nrows = 24264,stringsAsFactors = FALSE)
countGenes2 <- read.table(countGenespath, skip = 24264,stringsAsFactors = FALSE )

countGenes$Symbol <- countGenes[,1]
tmp <- match(countGenes$Symbol,allgenes_info$external_gene_id)
countGenes$geneID <- rownames(allgenes_info)[tmp]
countGenes <- countGenes[!is.na(countGenes$geneID),]
countGenes <- countGenes[!duplicated(countGenes$geneID),] 
rownames(countGenes) <- countGenes$geneID
countGenes <- countGenes[,-c(length(countGenes[1,])-1,length(countGenes[1,]))]
countGenes <- countGenes[,8:length(countGenes[1,])]
COL_NAME <- c('A01','B01','C01','D01','E01','F01','G01','H01','A02','B02','C02','D02','E02','F02','G02','H02',	
              'A03','B03','C03','D03','E03','F03','G03','H03','A04','B04','C04','D04','E04','F04','G04','H04','A05','B05',
              'C05','D05','E05','F05','G05','H05','A06','B06','C06','D06','E06','F06','G06','H06','A07','B07','C07','D07',
              'E07','F07','G07','H07','A08','B08','C08','D08','E08','F08','G08','H08','A09','B09','C09','D09','E09','F09',
              'G09','H09','A10','B10','C10','D10','E10','F10','G10','H10','A11','B11','C11','D11','E11','F11','G11','H11',	
              'A12','B12','C12','D12','E12','F12','G12','H12')
colnames(countGenes) <- COL_NAME
countGenes[,1:length(countGenes[1,])] <- lapply(countGenes[,1:length(countGenes[1,])],as.numeric)
well_quality_cell <- c("A03","B05","B06","B07","B09","B10","B12","C01","C02","C06","C08","A04",
                       "C09","C10","C11","C12","D01","D05","D06","D07","D11","E07","A05","E09",
                       "E12","F05","F10","G04","G05","G09","H04","H05","H06","A07","H08","H10",
                       "A08","A10","A11","B01","B04")
countGenes <- countGenes[,well_quality_cell]
colnames(countGenes) <- paste('SC_2i_',seq(1,41),sep = '')



# quality control
countGenes <- clean.counts(countGenes,min.detected = 2,min.reads = 1)
erccCount <- clean.counts(as.matrix(erccCount),min.detected = 2,min.reads = 1,min.lib.size = 20)
erccNumber <- erccNumber[rownames(erccCount),]
erccSeqID <- sort(rownames(erccNumber))
plot( log10(apply(erccCount, 1, var)), log10((apply(erccCount, 1, var)/apply(erccCount, 1, mean)^2)),xlim=c(-2,5),ylim=c(-2,2))
lines(seq(from=-2.0,to=5,length.out=100),-0.5*seq(from=-2.0,to=5,length.out=100),lwd=2)
# discard cells with fewer than 500 sequenced transcripts for ERCC spike-ins and 10,000 sequenced transcripts for endogenous genes.
removeCells = colSums(erccCount)>500 & colSums(countGenes)>10000
erccCount = erccCount[, removeCells]
countGenes = countGenes[, removeCells]

dim(countGenes)
head(countGenes[1:6, 1:6])


# Adjusting for batch effects
#cellCondition = "SC_serum"
cellCondition = "SC_2i"
nCount = selectCells(erccCount, countGenes, cellCondition, 1, 41, erccSeqID, erccNumber)
nCountSpikes = nCount[[1]]
numberSpikes = nCount[[2]]
nCountGenes = nCount[[3]]

sizeFactorMatrix = matrix(1, nrow=nrow(nCountSpikes), ncol=ncol(nCountSpikes))
gammaThetaEstimate = estimateGammaTheta(nCountSpikes, numberSpikes, sizeFactorMatrix)

EGammaThetaSC2i = gammaThetaEstimate$gammaTheta[[1]]
nCountGenesSC2i = nCountGenes / EGammaThetaSC2i
nCountSpikesSC2i = nCountSpikes / EGammaThetaSC2i



#Quantifying biological noise
sizeFactorMatrixSC2i = matrix(1, nrow=nrow(nCountSpikesSC2i),
                              ncol=ncol(nCountSpikesSC2i))
sizeFactorMatrixSC2iGenes = matrix(1, nrow=nrow(nCountGenesSC2i),ncol=ncol(nCountGenesSC2i))

noiseEstimateSC2i = estimateBiologicalVariance(nCountGenesSC2i, nCountSpikesSC2i,
                                               sizeFactorMatrixSC2i, numberSpikes,
                                               sizeFactorMatrixSC2iGenes)

# remove zero
noiseEstimateSC2i <- noiseEstimateSC2i[noiseEstimateSC2i$predictedCount>0.01,]
noiseEstimateSC2i <- noiseEstimateSC2i[noiseEstimateSC2i$predictedBVariance>0.01,]                                   
dim(noiseEstimateSC2i)

# plot some figure
plot(log10(noiseEstimateSC2i$predictedCount),log10(sqrt(noiseEstimateSC2i$predictedCV2)),
     xlim = c( -2, 3 ),ylim = c(-1.5, 1))
plot(log10(noiseEstimateSC2i$predictedCount),log10((noiseEstimateSC2i$predictedBVariance)),
     xlim = c( -2, 5 ),ylim = c(-4, 10))


gene_cellcycle_rank <- binx_nearesty(noiseEstimateSC2i$predictedCount, noiseEstimateSC2i$predictedBVariance, bin_size=101)
plot( log10(noiseEstimateSC2i$predictedCount), (gene_cellcycle_rank),xlim = c( -2, 5 ))

names(gene_cellcycle_rank) <- rownames(noiseEstimateSC2i)
gene_cellcycle_rank <- gene_cellcycle_rank[!is.na(gene_cellcycle_rank)]
length(gene_cellcycle_rank)
noiseEstimateSC2i <- noiseEstimateSC2i[names(gene_cellcycle_rank),]
mESCs_scRNA_removegene <- noiseEstimateSC2i


# Technical noise fit    
EVGammaThetaEstimate = estimateEVGammaTheta(nCountSpikesSC2i,
                                            numberSpikes, sizeFactorMatrixSC2i)
EGamma = EVGammaThetaEstimate$EGamma
ETheta = EVGammaThetaEstimate$ETheta
E2Gamma = EVGammaThetaEstimate$E2Gamma
E2Theta = EVGammaThetaEstimate$E2Theta
VGamma = EVGammaThetaEstimate$VGamma
VTheta = EVGammaThetaEstimate$VTheta


par(cex.axis=1, cex.lab=1)
plot( NULL, xaxt="n",
      log="xy", xlim = c( 1e-2, 1e5 ), ylim = c(0.01, 100),
      xlab = "Estimated number of transcripts per cell", ylab = "Squared CV" )
axis( 1, 10^(-2:5), c("0.01", "0.1", "1", "10", "100", "1000",
                      expression(10^4), expression(10^5)) )
points(rowMeans(nCountGenesSC2i) / (EGamma*ETheta),
       apply(nCountGenesSC2i, 1, var)/rowMeans(nCountGenesSC2i)^2,
       pch=20, cex=0.7, col="darkgray")
points(numberSpikes[,1], apply(nCountSpikesSC2i, 1, var)/rowMeans(nCountSpikesSC2i)^2,
       pch=20, cex=0.7, col="darkblue",lwd=2)
xg = 10^seq(-2, 5, length.out=1000)
predictedVariance = EGamma*ETheta*xg + (VGamma+E2Gamma)*(ETheta-(E2Theta+VTheta))*xg +(VGamma+E2Gamma)*VTheta*xg^2 + (E2Theta*VGamma)*xg^2
lines(xg, predictedVariance/(xg*EGamma*ETheta )^2, col="indianred", lwd=3)

save(mESCs_scRNA_removegene, gene_cellcycle_rank,
     file=paste(work_directory,'/data/mESCs_scRNA_islam_deconv', '.RData', sep=''))

