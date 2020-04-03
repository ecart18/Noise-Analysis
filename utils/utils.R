
binx_nearesty<-function(var_x, var_y, bin_size=101){
  temp_names<-names(var_y)  
  temp_order<-1:length(var_x)   
  temp_order<-temp_order[order(var_x)]  
  var_y<-var_y[order(var_x)]  
  var_x<-var_x[order(var_x)]  
  var_xi<-lapply(as.list(1:length(var_x)), function(x){c(x, seq(x-(bin_size-1)/2, x+(bin_size-1)/2)[-((bin_size-1)/2+1)])})
  var_ybr<-sapply(lapply(lapply(var_xi, function(x){var_y[x[x>0]]}), rank), '[', 1)
  var_ybr[sapply(var_xi, min)<1 | sapply(var_xi, max)>length(var_xi)]<-NA
  var_ybr<-var_ybr[order(temp_order)]
  names(var_ybr)<-temp_names
  var_ybr
}


detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}


clean.counts <- function(counts, min.lib.size = 1.8e3, min.reads = 10, min.detected = 5) {
  # filter out low-gene cells
  counts <- counts[, colSums(counts>0)>min.lib.size]
  # remove genes that don't have many reads
  counts <- counts[rowSums(counts)>min.reads, ]
  # remove genes that are not seen in a sufficient number of cells
  counts <- counts[rowSums(counts>0)>min.detected, ]
  return(counts)
}



# local regression
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family="binomial"), ...)
}

#Extract values from numeric input vector when index vector equals specified value
index_vect<-function(num_vect, index_vect, index_val){
  num_vect[which(index_vect==index_val)]
}


#Wrapper for Mann-Whitney U test
#Returns named vector with effect size (AUC) and p-value obtained from coin package
mann_whitney_U_wrapper_single<-function(vals1, vals2){
  if(length(vals1)==0 | length(vals2)==0){
    temp_test<-c(NA, NA)
  }else{
    g = factor(c(rep("GroupA", length(vals1)), rep("GroupB", length(vals2))))
    v = c(vals1, vals2)
    g_wt<-wilcox_test(v ~ g)
    temp_test<-c(wilcox.test(vals1, vals2)$statistic/(as.numeric(length(vals1))*as.numeric(length(vals2))),
                 pvalue(g_wt),length(vals1))
  }
  names(temp_test)<-c('effect_size', 'p_value','NumTarget')
  temp_test
}


#Wrapper for Mann-Whitney U test
#Returns named vector with effect size (AUC) and p-value obtained from coin package
mann_whitney_U_wrapper<-function(vals1, vals2){
  if(length(vals1)==0 | length(vals2)==0){
    temp_test<-c(NA, NA, NA)
  }else{
    g = factor(c(rep("GroupA", length(vals1)), rep("GroupB", length(vals2))))
    v = c(vals1, vals2)
    g_wt<-wilcox_test(v ~ g)
    
    # datatmp <- data.frame( c( rep(1,length(vals1)), rep(0, length(vals2))),
    #                        c(vals1,vals2))
    # colnames(datatmp) <- c('Group','Rank')
    # temp_test<-c(as.numeric(roc(datatmp$Group, datatmp$Rank,direction = "<")$auc),
    #              pvalue(g_wt),length(vals1))
    
    
    temp_test<-c(wilcox.test(vals1, vals2)$statistic/(as.numeric(length(vals1))*as.numeric(length(vals2))),
                 pvalue(g_wt),length(vals1))
    
  }
  names(temp_test)<-c('effect_size', 'p_value','NumTarget')
  temp_test
}


#Perform Mann-Whitney U test on numeric vectors in input lists (element indeces correspond)
#If background list not supplied, perform Mann-Whitney U test on each pair of numeric vectors in input list (first element in pair considered test vector)
#Returns data frame with effect size (AUC) and p-value obtained from coin package
mann_whitney_U_list<-function(test_list, background_list=NULL,sum_exp=NULL){
  if(is.null(background_list)){
    background_list<-test_list[seq(2, length(test_list), 2)]
    test_list<-test_list[seq(1, length(test_list), 2)]
  }
  temp_test<-as.data.frame(t(mapply(mann_whitney_U_wrapper, test_list, background_list)))
  temp_test$FDR<-p.adjust(temp_test$p_value,method = "BH")
  if(length(sum_exp) == dim(temp_test)[1]){temp_test$sum_exp <- sum_exp}
  temp_test
}

#Test biased values for all gene sets (supplied as list)
#gene_sets_list either contains binary (double) vector matching num_vect or character vector matching names of num_vect
genesets_biased<-function(gene_sets_list, num_vect){
  if(typeof(gene_sets_list[[1]])=='double'){
    temp_values<-as.list(as.data.frame(matrix(num_vect, nrow=length(num_vect), ncol=length(gene_sets_list), byrow=F)))  
    temp_in_list<-mapply(index_vect, temp_values, gene_sets_list, 1)  
    temp_out_list<-mapply(index_vect, temp_values, gene_sets_list, 0)  
    names(temp_in_list)<-names(gene_sets_list)   
    names(temp_out_list)<-names(gene_sets_list)
    mann_whitney_U_list(temp_in_list, temp_out_list)    
  }else if(typeof(gene_sets_list[[1]])=='character'){
    temp_in_list<-lapply(gene_sets_list, function(x){num_vect[names(num_vect) %in% x]})
    temp_out_list<-lapply(gene_sets_list, function(x){num_vect[!names(num_vect) %in% x]})
    mann_whitney_U_list(temp_in_list, temp_out_list)
  }
}



#Test biased values for all gene sets (supplied as list)
#gene_sets_list either contains binary (double) vector matching num_vect or character vector matching names of num_vect
genesets_biased_NonTargetCtrl<-function(gene_sets_list, num_vect, target_matrix, nontarget){
  if(typeof(gene_sets_list[[1]])=='double'){
    temp_values<-as.list(as.data.frame(matrix(num_vect, nrow=length(num_vect), ncol=length(gene_sets_list), byrow=F)))  
    temp_in_list<-mapply(index_vect, temp_values, gene_sets_list, 1)  
    
    num_vect_tmp <- num_vect[rownames(nontarget)[rownames(nontarget) %in% names(num_vect)]]
    temp_out_list <- as.list(as.data.frame(matrix(num_vect_tmp, nrow=length(num_vect_tmp), ncol=length(gene_sets_list), byrow=F)))
    
    
    names(temp_in_list)<-names(gene_sets_list)   
    names(temp_out_list)<-names(gene_sets_list)
    mann_whitney_U_list(temp_in_list, temp_out_list,sum_exp=NULL)    #  
  }else if(typeof(gene_sets_list[[1]])=='character'){
    temp_in_list<-lapply(gene_sets_list, function(x){num_vect[names(num_vect) %in% x]})
    num_vect_tmp <- num_vect[rownames(nontarget)[rownames(nontarget) %in% names(num_vect)]]
    temp_out_list <- as.list(as.data.frame(matrix(num_vect_tmp, nrow=length(num_vect_tmp), ncol=length(gene_sets_list), byrow=F)))
    mann_whitney_U_list(temp_in_list, temp_out_list)
  }
}
