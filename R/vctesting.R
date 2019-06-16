######################################################################
## Editor: Shiquan Sun and Jiaqiang Zhu
## Date: 2018-10-28 22:03:36
## Modified: 2019-4-26 11:24:27
## Affiliation: University of Michigan
######################################################################


#' Testing multiple kernel matrices
#' 
#' @param object SPARK object
#' @param kernel_mat A list to store the pre-defined kernel matrix, 
#' @param check_positive Check the kernel matrix is positive or not
#' @param verbose Output fitting information
#' @export
spark.test <- function(object, kernel_mat = NULL, check_positive = TRUE, verbose = TRUE) {
  
	if(!is.null(kernel_mat) & !is.list(kernel_mat)){ stop("SPARK.TEST::kernel_mat must be a list, please provide list type!")}
	# check one kernel at a time
    res_kernel <- list()
    res_pval <- NULL
	if(is.null(kernel_mat)){ #calculate kernel matrix using location information
		# Euclid distance, and compute the range of parameter l
		ED <- as.matrix(dist(object@location[ ,1:2]))
		lrang <- ComputeGaussianPL(ED, compute_distance=FALSE)[3:7]
    
		for(ikernel in c(1:5) ){
			# Gaussian kernel
			cat(paste0("## testing Gaussian kernel: ",ikernel,"...\n"))
			kernel_mat <- exp(-ED^2/(2*lrang[ikernel]^2))
			object <- spark.test_each(object, kernel_mat=kernel_mat, check_positive=check_positive, verbose=verbose)
			res_pval <- cbind(res_pval, object@res_stest$sw)
			res_kernel <- c(res_kernel, object@res_stest)
			rm(kernel_mat)

			# Periodic kernel
			cat(paste0("## testing Periodic kernel: ",ikernel,"...\n"))
			kernel_mat <- cos(2*pi*ED/lrang[ikernel])
			object <- spark.test_each(object, kernel_mat=kernel_mat, check_positive=check_positive, verbose=verbose)
			res_pval <- cbind(res_pval, object@res_stest$sw)
			res_kernel <- c(res_kernel, object@res_stest)
			rm(kernel_mat)
		}# end for ikernel
		colnames(res_pval) <- paste0(c("GSP","COS"), rep(1:5,each=2))
		rownames(res_pval) <- rownames(object@counts)
	}else{# kernel_mat is a list provided by user
		for(ikernel in 1:length(kernel_mat) ){
			# pre-defined kernels
			object <- spark.test_each(object, kernel_mat=kernel_mat[[ikernel]], check_positive=check_positive, verbose=verbose)
			res_pval <- cbind(res_pval, object@res_stest$sw)
			res_kernel <- c(res_kernel, object@res_stest)
		}# end for ikernel
		colnames(res_pval) <- paste0("kernel", 1:length(kernel_mat) )
		rownames(res_pval) <- rownames(object@counts)
	}# end fi
    
    object@res_stest <- res_kernel
    ## summarize ten pvalues into one
    combined_pvalue <- CombinePValues(res_pval)
    object@res_mtest <- data.frame(res_pval, combined_pvalue = combined_pvalue,  adjusted_pvalue = p.adjust(combined_pvalue, method="BY") )
	# return the results
    return(object)
}# end function score test



#' Testing one kernel matrix to identify spatial pattern
#' 
#' @param object SPARK object
#' @param kernel_mat The kernel matrix
#' @param check_positive Check the kernel matrix is positive or not
#' @param verbose Output fitting information
#' @export
#' 
spark.test_each <- function(object, kernel_mat, check_positive = FALSE, verbose = TRUE) {
    
  # The number of core used in testing step
  num_core <- object@num_core
  # using parallel to correct
  #library(doSNOW)
  if(num_core > 1){
       if(num_core > detectCores()){warning("SPARK:: the number of cores you're setting is larger than detected cores!");num_core = detectCores()}
  }#end fi
  cl <- makeCluster(num_core)
  registerDoSNOW(cl)
	
  num_cell <- ncol(object@counts)
  if(check_positive){# kernel matrix should be positive definition matrix
    # need to check positive definition before tesing
	# if(verbose){cat(paste("SPARK.TEST:: Checking the kernel matrix. \n"))}
    # eig <- eigen(kernel_mat)
    # eigval <- eig$value
    # eigvector <- eig$vectors
    # if(any(eigval<1e-8)){ 
    ### warning("SPARK.TEST::the kernel matrix is singular, it has been modified!")
      # eigval[eigval<1e-8] <- 1e-8
      # kernel_mat <- eigvector%*%diag(eigval)%*%t(eigvector)
    # }# end fi
    # rm(eig)
    # rm(eigval)
    # rm(eigvector)
	res <- SysMatEigen(kernel_mat)
	kernel_mat <- res$kernel_mat
	# if(num_cell < 20000){
		# res <- SysMatEigen(kernel_mat)
		# kernel_mat <- res$kernel_mat
	# }else{
		# kernel_mat[which(kernel_mat<1e-6)] <- 0
		# kernel_mat <- as(kernel_mat, "dgCMatrix")
		# res <- SparseSysMatEigen(kernel_mat, round(0.04*num_cell) ) # extract 4% PCs
		# kernel_mat <- res$kernel_mat
	# }# end fi
	rm(res)
  }# end fi check
  
  # number of genes to testing
  num_gene_test <- length(object@res_vc)
  #========================================
  # using parallel to test variance components
  #registerDoParallel(cores = num_core)
  pb <- txtProgressBar(max = num_gene_test, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
	
  res_test <-foreach(ig = 1:num_gene_test, .combine=rbind, .options.snow=opts)%dopar%{
    if(verbose) {cat(paste("NO. Gene = ",ig,"\n"))}
    model1 <- object@res_vc[[ig]]
    
    davies_sw <-NA
    converged <- FALSE
    if((class(model1) != "try-error") & (!any(is.na(model1$Y)))){
      converged <- model1$converged
      # extract the analyzed cells
      cov_mat <- kernel_mat[model1$analy_indx, model1$analy_indx]
      if(ncol(model1$X)>1){
        rest <- ComputeTestQuantRcpp_cov(model1$Y, model1$Py, model1$X, cov_mat, model1$D^2, model1$theta)
      }else{
        rest <- ComputeTestQuantRcpp_nocov(model1$Y, model1$Py, cov_mat, model1$D^2, model1$theta)
      }# end fi
      
      #scaledchisq.pvalue 	<- pchisq(rest$S0/rest$kk, rest$df, lower.tail = F)
      ##--------------------
      ## davies_sw
      newInfoM 	<- rest$infoMp1
      # calculate the scale parameter and degree of freedom
      ee_sw <- rest$ee
      kk_sw <- newInfoM/(2*ee_sw)
      df_sw <- 2*ee_sw^2/(newInfoM)
      davies_sw	<- pchisq(rest$S0/kk_sw, df_sw, lower.tail = F)
      if(verbose) {cat(paste("SPARK.SCORE::SW pvalue 1 = ", davies_sw,"\n"))}
    } # end fi
    # to return
    #return( data.frame(geneid = names(object@res_vc)[ig], p_score = scaledchisq.pvalue, p_davies = p_davies, sw=davies_sw, converged=converged) )
    return( data.frame(geneid = names(object@res_vc)[ig], sw=davies_sw, converged=converged) )
  }# end parallel foreach
  close(pb)
  stopCluster(cl)
  #######
  object@res_stest <- res_test
  rm(res_test)
  return(object)
}# end function score test for each kernel

#########################################
#             CODE END                  #
#########################################




