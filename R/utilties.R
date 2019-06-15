########################################
# Edited by Shiquan Sun
# Date: 2018-12-29 07:48:59
# Modified: 2019-4-28 11:53:55
##########################################

#' Compute Gaussian kernel matrix
#' @param data1 The data set
#' @param data2 The second data set
#' @param gamma The parameter for Guassian kernel
#' @export
# the suggested gamma is taking 0.01, 0.05, 0.1, 0.2, 0.5
ComputeGaussianL <- function(data1, data2=NULL, gamma=0.1){
	data1 <- as.matrix(data1)
	num_cell <- nrow(data1)
	if(is.null(data2)){
		data2 <- data1
	}# end fi
	
	dist_mat <- matrix(0, ncol = num_cell, nrow = num_cell)
	kernel_para <- rep(0, nrow(data1))
	for(k in 1:nrow( data1)){ #
		# dist_mat[k,]  <-  rowSums( (data2 -  kronecker(matrix(rep(1,nrow(data2)),ncol=1), matrix(data1[k,],nrow = 1), FUN = "*")  )^2 )
		dist_mat[k,] <- rowSums( abs(data2 - kronecker(matrix(rep(1,nrow(data2)),ncol=1), matrix(data1[k,],nrow = 1), FUN = "*")  ) )
		PP <- -dist_mat[k,]/log(gamma)
		#kernel_para[k] <- min( PP[which(abs(PP)>0)] )
		kernel_para[k] <- mean( PP[which(abs(PP)>0)] )
	}# end for
	
	# return cell length parameter
	return(mean(kernel_para))
}# end func


#' Filtering out the genes that are lowly expressed acorss cells or the cells that are a few number of genes expressed
#' @param counts Raw gene expression counts, p x n matrix
#' @param prectage_cell The prectage of cells that are expressed
#' @param min_total_counts Minimum counts for each cell
#' @export
FilteringGenesCells <- function(counts, prectage_cell=0.1, min_total_counts=10){
  idx <- which( apply(counts ,1 ,function(x){sum(x!=0)})>=(floor(prectage_cell*ncol(counts))) )
  # filtering out genes
  counts <- counts[idx,]
  # filtering out cells
  counts <- counts[ ,which(apply(counts,2,sum)>min_total_counts)]
  return(counts)
}# end func




## Calculate Gaussian Parameter L based on Eculidean Distance
#' @param X Cell corrdinates matrix n x 2 or kernel matrix computed already
#' @param distance Compute the distance matrix using generic function dist
#' @export
ComputeGaussianPL <- function(X, compute_distance=TRUE){
  if(compute_distance){
    if(ncol(X)<2){stop("X has to be a coordinate matrix with number of column greater than 1")}
    D <- as.matrix(dist(X))  
  }else{
    D <- X
  }# end fi
  
  Dval <- unique(as.vector(D))
  Dval.nz <- Dval[Dval>1e-8]
  lmin <- min(Dval.nz)/2
  lmax <- max(Dval.nz)*2
  lrang <- 10^(seq(log10(lmin),log10(lmax),length.out=10))
  return(lrang)
}# end func


## Summarized the multiple p-values via Cauchy combination rule
#' @param pvalues Pvalues from multiple kernels
#' @export
CombinePValues <- function(pvalues){
  Cstat <- tan((0.5-pvalues)*pi)
  Cbar <- apply(Cstat, 1, mean)
  combined_pval <- 1.0/2.0 - atan(Cbar)/pi
  combined_pval[which(combined_pval<=0)] <- 5.55e-17
  return(combined_pval)
}# end func



