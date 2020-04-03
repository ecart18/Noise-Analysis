##### Internal R Functions

# select cells of a given cell condition (cellCondition) (a <= cell index <= b)
selectCells <- function(erccCount, countGenes, cellCondition, a, b, erccSeqID, erccNumber) {
  splitColNames = strsplit(colnames(erccCount), split="_")
  indexColumns = sapply(1:ncol(erccCount), function(x) {
    if (paste(splitColNames[[x]][1], splitColNames[[x]][2], sep="_")==cellCondition &
          as.numeric(splitColNames[[x]][3]) >= a &
          as.numeric(splitColNames[[x]][3]) <= b) {
      TRUE
    }
    else {
      FALSE
    }
  } )
  nCountSpikes = matrix(0, length(erccSeqID), ncol(erccCount[,indexColumns]))
  rownames(nCountSpikes) = as.character(erccSeqID)
  nCountSpikes[rownames(erccCount),] = as.matrix(erccCount[,indexColumns])
  numberSpikes = erccNumber[rownames(nCountSpikes),2]
  numberSpikes = cbind(numberSpikes, numberSpikes)
  
  splitColNames = strsplit(colnames(countGenes), split="_")
  indexColumns = sapply(1:ncol(countGenes), function(x) {
    if (paste(splitColNames[[x]][1], splitColNames[[x]][2], sep="_")==cellCondition &
          as.numeric(splitColNames[[x]][3]) >= a &
          as.numeric(splitColNames[[x]][3]) <= b) {
      TRUE
    }
    else {
      FALSE
    }
  } )
  nCountGenes <- as.matrix(countGenes[,indexColumns])
  list(nCountSpikes, numberSpikes, nCountGenes)
}

# Estimate gamma and theta
estimateGammaTheta <- function(nCountSpikes, numberSpikes, sizeFactorMatrix) {
  fitData = data.frame(Xi=numberSpikes[,1], Ki=rowMeans(nCountSpikes))
  fit1 <- lm(Ki ~ 0 + Xi, data=fitData) 
  gammaTheta = coefficients(fit1)    
  PropCell = sapply(1:nrow(nCountSpikes), function(x) {
    sum(nCountSpikes[x,]>0)/ncol(nCountSpikes)
  })
  fitData = data.frame(Xi=numberSpikes[,1], Y=PropCell, Ai=rowMeans(sizeFactorMatrix))
  initialTheta = mean(fitData$Y[fitData$Xi>0.5 & fitData$Xi<5])
  if (initialTheta==0 | is.na(initialTheta)) {
    initialTheta = 0.01
  }    
  fit2 <- nlsLM(Y ~ 1 - (1 - Theta + Theta*exp(-(gammaTheta/Theta)*Ai))^Xi, data=fitData, start=list(Theta=initialTheta),
                lower=0, upper=1, control=nls.control(warnOnly=TRUE))
  
  
  Theta = coefficients(fit2)
  list(gammaTheta=gammaTheta, Theta=Theta)
}

# Initializing V[gamma] and V[theta]
estimateVGammaThetaInitial <- function(nCountSpikes, numberSpikes, sizeFactorMatrix) {
  old.o = options()
  options(warn=FALSE)
  gammaThetaSample = aaply(1:ncol(nCountSpikes), 1, function(x) {
    fitData = data.frame(Xi=numberSpikes[,1], Ki=nCountSpikes[,x])
    fit1 <- lm(Ki ~ 0 + Xi, data=fitData) 
    gammaTheta = coefficients(fit1)
    
    fitData = data.frame(Xi=numberSpikes[,1], Y=nCountSpikes[,x]>0, Ai=sizeFactorMatrix[,x])
    
    fit2 <- nlsLM(Y ~ 1 - (1 - Theta + Theta*exp(-(gammaTheta/Theta)*Ai))^Xi, data=fitData,
                  start=list(Theta=0.4), lower=0, upper=1, control=nls.control(warnOnly=TRUE))    
    c(coefficients(fit2), gammaTheta/coefficients(fit2))}, .expand=FALSE)
  gammaThetaSample[gammaThetaSample[,1]==1,1] = 0.9999
  
  abTheta = try(fitdistr(gammaThetaSample[,1], "beta", list(shape1=0.5, shape2=0.5)), silent=TRUE)
  if (class(abTheta) == "try-error") {
    VTheta = var(gammaThetaSample[,1], na.rm=TRUE)
  } else {
    aTheta = abTheta$estimate[[1]]
    bTheta = abTheta$estimate[[2]]
    VTheta = (aTheta*bTheta) / ( (aTheta+bTheta)^2*(aTheta+bTheta+1))
  }
  abGamma = try(fitdistr(gammaThetaSample[,2], "gamma"), silent=TRUE)
  if (class(abGamma) == "try-error") {
    VGamma = var(gammaThetaSample[,2], na.rm=TRUE)
  } else {
    aGamma = abGamma$estimate[[1]]
    bGamma = abGamma$estimate[[2]]
    VGamma = aGamma/bGamma^2
  }
  options(old.o)
  c(VGamma, VTheta)
}

# Estimate E[gamma], E[theta], V[Gamma] and V[Theta]
estimateEVGammaTheta <- function(nCountSpikes, numberSpikes, sizeFactorMatrix) {
  gammaThetaEstimate = estimateGammaTheta(nCountSpikes, numberSpikes, sizeFactorMatrix)
  EGamma = gammaThetaEstimate$gammaTheta[[1]] / gammaThetaEstimate$Theta[[1]] 
  ETheta = gammaThetaEstimate$Theta[[1]] 
  E2Gamma = EGamma^2
  E2Theta = ETheta^2
  
  varianceGammaThetaEstimate = estimateVGammaThetaInitial(nCountSpikes, numberSpikes, sizeFactorMatrix)
  VGamma = varianceGammaThetaEstimate[1]
  VTheta = varianceGammaThetaEstimate[2]  
  
  fitData = data.frame(Xi=numberSpikes[,1], Y=apply(nCountSpikes,1,var)/rowMeans(nCountSpikes), Bi=rowMeans(1/sizeFactorMatrix))
  
  fit <- nlsLM(Y ~ Bi + (VGamma+E2Gamma)/EGamma*(1-(E2Theta+VTheta)/ETheta) + 
                 (VGamma+E2Gamma)*VTheta*Xi/(EGamma*ETheta) + (ETheta*VGamma)*Xi/EGamma, data=fitData,
               start=list(VGamma=VGamma, VTheta=VTheta), lower=c(0,0))  
  VGamma=coefficients(fit)[[1]]
  VTheta=coefficients(fit)[[2]]
  # second optimization for robust estimates
  fit <- nlsLM(Y ~ Bi + (VGamma+E2Gamma)/EGamma*(1-(E2Theta+VTheta)/ETheta) + 
                 (VGamma+E2Gamma)*VTheta*Xi/(EGamma*ETheta) + (ETheta*VGamma)*Xi/EGamma, data=fitData,
               start=list(VGamma=VGamma, VTheta=VTheta), lower=c(0,0))  
  VGamma=coefficients(fit)[[1]]
  VTheta=coefficients(fit)[[2]]
  list(EGamma=EGamma, ETheta=ETheta, E2Gamma=E2Gamma, E2Theta=E2Theta, VGamma=VGamma, VTheta=VTheta)
}

# Estimate biological variance
estimateBiologicalVariance <- function(nCountGenes, nCountSpikes, sizeFactorMatrix, numberSpikes, sizeFactorMatrixGenes) {
  EVGammaThetaEstimate = estimateEVGammaTheta(nCountSpikes, numberSpikes, sizeFactorMatrix)
  EGamma = EVGammaThetaEstimate$EGamma
  ETheta = EVGammaThetaEstimate$ETheta
  E2Gamma = EVGammaThetaEstimate$E2Gamma
  E2Theta = EVGammaThetaEstimate$E2Theta
  VGamma = EVGammaThetaEstimate$VGamma
  VTheta = EVGammaThetaEstimate$VTheta
  
  Bi = rowMeans(1/sizeFactorMatrixGenes)  
  predictedCount = rowMeans(nCountGenes) / (EGamma*ETheta)  
  totalVariance = apply(nCountGenes, 1, var)
  predictedBVariance = totalVariance - (Bi*EGamma*ETheta*predictedCount + (VGamma+E2Gamma)*(ETheta-(E2Theta+VTheta))*predictedCount + 
                                          (VGamma+E2Gamma)*VTheta*predictedCount^2 + (E2Theta*VGamma)*predictedCount^2)
  fracBVariance = predictedBVariance / totalVariance
  fracBVariance[fracBVariance<0] = 0
  predictedBVariance = predictedBVariance / ( (VGamma+E2Gamma)*(VTheta+E2Theta) )
  predictedBVariance[predictedBVariance<0] = 0
  predictedCV2 = predictedBVariance/predictedCount^2
  predictedFano = predictedBVariance/predictedCount
  data.frame(predictedCount=predictedCount, predictedBVariance=predictedBVariance,
             predictedCV2=predictedCV2, predictedFano=predictedFano, fracBVariance=fracBVariance, row.names=names(predictedCount))
}

repmat <- function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx <- dim(X)[1]
  nx <- dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}

# Estimate the effective length of SNP (See supplementary Note 7)
estimateSNPLength <- function(SNPAnnot, fragmentLength, readLength) {
  effectiveSize=apply(SNPAnnot, 1, function(snp) {
    trLength = as.numeric(snp[6])
    relPosition = as.numeric(snp[9])
    relPosition = min(relPosition, trLength - relPosition + 1)
    if (trLength < fragmentLength) {
      effectiveSize = 1
    } else if (relPosition <= fragmentLength & relPosition <= readLength) {
      effectiveSize = relPosition
    } else if (relPosition <= fragmentLength & relPosition > readLength) {
      effectiveSize = readLength + max(relPosition - fragmentLength + readLength, 0)
    } else {
      effectiveSize = 2*readLength
    }
    effectiveSize
  } 
  )
}

# Simulating single-cell data (See supplementary Note 6)
simulateCountEfficient <- function(Xi, ETheta, VTheta, EGamma, VGamma, Aij, nCell, BV, nSample) {
  if (VGamma > 0) {
    gammaShape = EGamma^2 / VGamma
    gammaScale = VGamma / EGamma
  }
  if (VTheta > 0) {
    thetaA = ETheta*( (ETheta*(1-ETheta))/VTheta - 1)
    thetaB = (1-ETheta)*( (ETheta*(1-ETheta))/VTheta - 1)
  }
  if (VGamma > 0) {
    gammaj = rgamma(nCell*nSample, shape=gammaShape, scale=gammaScale)
  } else {
    gammaj = rep(EGamma, nCell*nSample)
  }
  if (VTheta > 0) {
    thetaj = rbeta(nCell*nSample, shape1=thetaA, shape2=thetaB)
  } else {
    thetaj = rep(ETheta, nCell*nSample)
  }
  if (BV > 0) {
    Xi = rgamma(nCell*nSample, shape=rep(Xi^2/BV, nCell*nSample), scale=rep(BV/Xi, nCell*nSample))
  } else {
    Xi = rep(Xi, nCell*nSample)
  }   
  Zij = rbinom(nCell*nSample, round(Xi), thetaj)
  Kij = rpois(nCell*nSample, gammaj*rep(Aij, nSample)*Zij)  
  Kij = matrix(Kij, nSample, nCell, byrow=TRUE)
}

simulateCountGenes <- function(Xi, ETheta, VTheta, EGamma, VGamma, Aij, nCell, BV) {
  nGenes = length(Xi)
  if (VGamma > 0) {
    gammaShape = EGamma^2 / VGamma
    gammaScale = VGamma / EGamma
  }
  if (VTheta > 0) {
    thetaA = ETheta*( (ETheta*(1-ETheta))/VTheta - 1)
    thetaB = (1-ETheta)*( (ETheta*(1-ETheta))/VTheta - 1)
  }
  if (VGamma > 0) {
    gammaj = rgamma(nCell*nGenes, shape=gammaShape, scale=gammaScale)
  } else {
    gammaj = rep(EGamma, nCell*nGenes)
  }
  if (VTheta > 0) {
    thetaj = rbeta(nCell*nGenes, shape1=thetaA, shape2=thetaB)
  } else {
    thetaj = rep(ETheta, nCell*nGenes)
  }
  if (sum(BV) > 0) {
    Yi = rgamma(nCell*nGenes, shape=Xi^2/BV, scale=BV/Xi)
  } else {
    Yi = rep(Xi, nCell)
  }    
  Zij = rbinom(nCell*nGenes, round(Yi), thetaj)
  Kij = rpois(nCell*nGenes, gammaj*Aij*Zij)  
  Kij = matrix(Kij, nGenes, nCell, byrow=F)
  Kij[Xi==0,] = 0
  Kij
}


readGenesSpikesCountTable <- function(fileName, fileDirectory, replicate, mtGenes, badIndex) 
{
  setwd(fileDirectory)
  cell = read.table(file=paste(fileName, "_", replicate, ".txt", sep=""), 
                    sep="\t", stringsAsFactor=FALSE, header=TRUE, check.names=FALSE)
  countCell = as.matrix(cell[,-1])
  rownames(countCell) = cell[,1]
  colnames(countCell) = sapply(colnames(countCell), function (x) unlist(strsplit(tail(unlist(strsplit(x, "_")),1), "\\."))[[1]])
  
  geneCount = colSums(countCell[1:38561,])
  nogeneCount = countCell[38654,]
  nomappedCount = countCell[38657,]
  exonProp = rbind(geneCount, nogeneCount, nomappedCount)
  exonProp = t(t(exonProp) / colSums(exonProp))
  mtCount = colSums(countCell[mtGenes[,1],])
  mtProp = mtCount / geneCount
  
  badIndex = union(badIndex, names(which(geneCount<500000 | exonProp[1,] < 0.5 | mtProp > 0.1))) 
  
  if (is.null(badIndex) || length(badIndex)==0) {
    countCell = countCell[1:38561,]
  } else {
    countCell = countCell[1:38561, !colnames(countCell) %in% as.character(badIndex)]
  }
}

qcCountTable <- function(count_cell, bad_index, mt_genes)
{
  gene_count = colSums(count_cell[1:38561,])
  nogene_count = count_cell[38654,]
  nomapped_count = count_cell[38657,]
  exon_prop = rbind(gene_count, nogene_count, nomapped_count)
  exon_prop = t(t(exon_prop) / colSums(exon_prop))
  mt_counts = colSums(count_cell[mt_genes[,1],])
  mt_prop = mt_counts / gene_count  
  bad_index_final = union(bad_index, which(gene_count<500000 | exon_prop[1,] < 0.5 | mt_prop > 0.1))
  count_cell = count_cell[1:38561, !colnames(count_cell) %in% as.character(bad_index_final)]
}

estimateBiologicalCVBootstrap <- function(nCountGenes, nCountSpikes, numberSpikes, genes, nB=100) {
  nCell = ncol(nCountSpikes)
  sizeFactorMatrix = matrix(1, nrow=nrow(nCountSpikes), ncol=nCell)
  EVGammaThetaEstimate = estimateEVGammaTheta(nCountSpikes, numberSpikes, sizeFactorMatrix)
  EGamma = EVGammaThetaEstimate$EGamma
  ETheta = EVGammaThetaEstimate$ETheta
  E2Gamma = EVGammaThetaEstimate$E2Gamma
  E2Theta = EVGammaThetaEstimate$E2Theta
  VGamma = EVGammaThetaEstimate$VGamma
  VTheta = EVGammaThetaEstimate$VTheta
  
  Bi = 1 # rowMeans(1/sizeFactorMatrix)
  predictedCVSE = sapply(genes, function(x) {
    predictedCV=sapply(1:nB, function(ii) {
      y = sample(nCountGenes[x,], nCell, replace=TRUE)
      if (mean(y) > 0) {
        predictedCount = mean(y) / (EGamma*ETheta)
        predictedBVariance = var(y) - (Bi*EGamma*ETheta*predictedCount + (VGamma+E2Gamma)*(ETheta-(E2Theta+VTheta))*predictedCount + 
                                         (VGamma+E2Gamma)*VTheta*predictedCount^2 + (E2Theta*VGamma)*predictedCount^2)
        predictedBVariance[predictedBVariance<0] = 0
        predictedBVariance = predictedBVariance / ( (VGamma+E2Gamma)*(VTheta+E2Theta) )
        predictedCV = sqrt(predictedBVariance/predictedCount^2)
      } else {
        predictedCV = 0
      }
    })
    predictedCV = predictedCV[predictedCV>0]
    c(mean(predictedCV), sd(predictedCV))
  })
}