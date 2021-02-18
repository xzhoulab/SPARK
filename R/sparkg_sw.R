#' @title Fast Test for the Gaussian version of SPARK through Satterthwaite Approximation
#' @description Efficient Test with Satterthwaite Approximation in Gaussian Models.
#' @param counts: A n x p gene expression matrix (sparseMatrix)
#' @param location: A n x d location matrix 
#' @param verbose: A boolean value whether to print details, for debug purpose
#' @return A list of estimated parameters
#' \item{stats}{A n x 10 matrix of test statistics for all kernels}
#' \item{res_stest}{A n x 10 matrix of P values for all kernels}
#' \item{res_mtest}{A n x 2 matrix of combined P values and BY-adjusted P values}
#' @export
sparkg_sw <- function(counts,location,verbose=FALSE){
  gene_rm_idx   <- which(sp_sums_Rcpp(counts,rowSums=TRUE)==0)
  if(length(gene_rm_idx)!=0){
    counts     <- counts[-gene_rm_idx,]
  }

  cell_rm_idx   <- which(sp_sums_Rcpp(counts)==0)
  if(length(cell_rm_idx)!=0){
    count_mat     <- counts[,-cell_rm_idx]
    location_mat  <- location[-cell_rm_idx,]
  }else{
    count_mat     <- counts
    location_mat  <- location
  }

  numGene     <- nrow(count_mat)
  numCell     <- ncol(count_mat)

  cat(paste("## ===== SPARK-G INPUT INFORMATION ==== \n"))
  cat(paste("## number of total samples:", numCell,"\n"))
  cat(paste("## number of total genes:", numGene,"\n"))

  scaled_count  <- NormalizeVSTfast(count_mat)
  ED            <- ED_cpp(as.matrix(location_mat))
  lrang         <- ComputeGaussianPL(ED,compute_distance=FALSE)

  res_sparkg <- list()
  for(iker in 1:10){
    cat("## Testing Kernel",iker,"\n")
    kernel_mat         <- exp(-ED^2/(2*lrang[iker]^2))
    res_sparkg[[iker]] <- sparkg_sw_sk(scaled_count,kernel_mat,verbose=verbose)
    rm(kernel_mat)
  }

  all_stats <- as.data.frame(sapply(res_sparkg,function(x){x[,1]}))
  all_pvals <- as.data.frame(sapply(res_sparkg,function(x){x[,2]}))

  rownames(all_pvals) <- rownames(all_stats) <- rownames(counts)
  colnames(all_pvals) <- colnames(all_stats) <- paste0("Gaussian",1:10)

  comb_pval   <- apply(all_pvals,1,ACAT)
  pBY         <- p.adjust(comb_pval,method="BY")
  joint_pval  <- cbind.data.frame(combinedPval=comb_pval,adjustedPval=pBY)

  final_res   <- list(stats=all_stats,res_stest=all_pvals,res_mtest=joint_pval)
  return(final_res)
}


#' @title Test for single kernel in the fast version of SPARK-G
#' @param Y: A n x p VST transformed gene expression matrix 
#' @param Kmat: A n x n kernel matrix of interest
#' @param verbose: A boolean value whether to print details, for debug purpose
#' @return A dataframe of test statistics and P values
#' @export
sparkg_sw_sk <- function(Y,Kmat,verbose=FALSE){
  num_gene    <- nrow(Y)
  num_cell    <- ncol(Y)
  covariates  <- as.matrix(rep(1,num_cell))
  Xdagger     <- pinv(covariates)
  SK          <- Kmat - covariates%*%(Xdagger%*%Kmat)

  ES          <- 0.5*sum(diag(SK))
  VS          <- 0.5*sum(t(SK)*SK) ## tr(SKSK)

  kk          <- VS/(2*ES)
  vv          <- 2*ES^2/VS

  sw_stat <- pval_vec <- c()
  for(igene in 1:num_gene){
    yy      <- Y[igene,]
    meanY   <- mean(yy)
    varY    <- var(yy)
    SY      <- yy - meanY
    S0      <- 0.5*crossprod(SY,crossprod(Kmat,SY))/varY

    sw_pval   <- pchisq(S0/kk, vv, lower.tail = F)
    pval_vec  <- c(pval_vec,sw_pval)
    sw_stat   <- c(sw_stat,S0)
    if(verbose){cat("Gene",igene," SW pvalue = ",sw_pval,"\n")}
  }
  sparkg_sk_out <- cbind.data.frame(score=sw_stat,pval=pval_vec)
  rownames(sparkg_sk_out) <- rownames(Y)
  return(sparkg_sk_out)
}


