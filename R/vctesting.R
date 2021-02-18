######################################################################
# Package: SPARK
# Version: 1.0.1
# Date   : 2018-10-23
# Modified: 2019-5-20 08:20:20;2019-8-10 19:24:50
# Title  : Count-based spatial model for identifying spatially variable genes
# Authors: S.Q. Sun, J.Q. Zhu, and X. Zhou
# Contacts: shiquans@umich.edu and jiaqiang@umich.edu 
#          University of Michigan, Department of Biostatistics
######################################################################


#' Testing multiple kernel matrices
#' 
#' @param object SPARK object
#' @param kernel_mat A list to store the pre-defined kernel matrix, 
#' @param check_positive Check the kernel matrix is positive or not
#' @param verbose Output fitting information
#' @export
spark.test <- function(object, 
						kernel_mat = NULL, 
						check_positive = TRUE, 
						verbose = TRUE) {
	
	## the variables input
	if(!is.null(kernel_mat) & !is.list(kernel_mat)){ stop("SPARK.TEST::kernel_mat must be a list, please provide list type!")}
	## testing one kernel at a time
    res_kernel <- list()
    res_pval <- NULL
	if(is.null(kernel_mat)){ #calculate kernel matrix using location information
		# Euclid distance, and compute the range of parameter l
		ED <- as.matrix(dist(object@location[ ,1:2]))
		lrang <- ComputeGaussianPL(ED, compute_distance=FALSE)[3:7]
		## total 10 kernels, i.e., 5 gaussian and 5 periodic kernels
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
			cat(paste0("## testing pre-defined kernel: ",ikernel,"...\n"))
			object <- spark.test_each(object, kernel_mat=kernel_mat[[ikernel]], check_positive=check_positive, verbose=verbose)
			res_pval <- cbind(res_pval, object@res_stest$sw)
			res_kernel <- c(res_kernel, object@res_stest)
		}# end for ikernel
		colnames(res_pval) <- paste0("kernel", 1:length(kernel_mat) )
		rownames(res_pval) <- rownames(object@counts)
	}# end fi
    
    object@res_stest <- res_kernel
    ## integrate ten p-values into one
	num_pval <- ncol(res_pval)
	num_gene <- nrow(res_pval)
	## compute the weight to p-values integrate
	if(is.null(object@weights)){
		weights <- matrix(rep(1.0/num_pval, num_pval*num_gene), ncol=num_pval )
	}else if(!is.matrix(object@weights)){
		weights <- as.matrix(object@weights)
	}else {
		weights <- object@weights
	}# end fi
    combined_pvalue <- CombinePValues(res_pval, weights=weights)
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
	## The number of core used in testing step
	num_core <- object@num_core
	## using parallel to correct
	#library(doSNOW)
	if(num_core > 1){
       if(num_core > detectCores()){warning("SPARK:: the number of cores you're setting is larger than detected cores!");num_core = detectCores()}
	}#end fi
	#cl <- makeCluster(num_core)
	#registerDoSNOW(cl)
	
	num_cell <- ncol(object@counts)
	
	## number of genes to testing
	num_gene_test <- length(object@res_vc)
	fit.model <- object@fit_model
	#========================================
	
	if(fit.model=="poisson"){
		if(check_positive){## kernel matrix should be positive definition matrix
			res <- SysMatEigen(kernel_mat)
			kernel_mat <- res$kernel_mat
			rm(res)
		}# end fi check
		## using parallel to test variance components
		registerDoParallel(cores = num_core)
		#pb <- txtProgressBar(max = num_gene_test, style = 3)
		#progress <- function(n) setTxtProgressBar(pb, n)
		#opts <- list(progress = progress)
		#res_test <-foreach(ig = 1:num_gene_test, .combine=rbind, .options.snow=opts)%dopar%{
		ig <- 1
		res_test <- foreach(ig = 1:num_gene_test, .combine=rbind)%dopar%{
		#res_test <- NULL
		#for(ig in 1:num_gene_test){
			if(verbose) {cat(paste("NO. Gene = ",ig,"\n"))}
			model1 <- object@res_vc[[ig]]
    
			davies_sw <- NA
			converged <- FALSE
			if((class(model1) != "try-error") & (!any(is.na(model1$Y)))){
				converged <- model1$converged
				## extract the analyzed cells
				cov_mat <- kernel_mat[model1$analy_indx, model1$analy_indx]
				if(ncol(model1$X)>1){
					rest <- ComputeTestQuantRcpp_cov(model1$Y, model1$Py, model1$X, cov_mat, model1$D^2, model1$theta)
				}else{
					rest <- ComputeTestQuantRcpp_nocov(model1$Y, model1$Py, cov_mat, model1$D^2, model1$theta)
				}# end fi
      
				#scaledchisq.pvalue 	<- pchisq(rest$S0/rest$kk, rest$df, lower.tail = F)
				##--------------------
				## davies_sw
				newInfoM <- rest$infoMp1
				## calculate the scale parameter and degree of freedom
				ee_sw <- rest$ee
				kk_sw <- newInfoM/(2*ee_sw)
				df_sw <- 2*ee_sw^2/(newInfoM)
				davies_sw <- pchisq(rest$S0/kk_sw, df_sw, lower.tail = F)
				if(verbose) {cat(paste("SPARK.SCORE::SW pvalue 1 = ", davies_sw,"\n"))}	
			} # end fi
			## to return
			res_each <- data.frame(geneid = names(object@res_vc)[ig], sw=davies_sw, converged=converged)
			
		}# end parallel foreach
	}else if(fit.model=="gaussian"){
		## required packages
		#library(pracma)
		#library(CompQuadForm)
		#check_positive <- FALSE
		if(check_positive){## kernel matrix should be positive definition matrix
			res <- SysMatEigen(kernel_mat)
			kernel_mat <- res$kernel_mat
			rm(res)
		}# end fi check
		zeros_threshold <- 10e-5
		norm_counts <- t(object@scaled_counts)
		covariates <- object@res_vc$covariates
		## check data format
		# if(class(kernel_mat) != "matrix"){
		if(!is(kernel_mat, "matrix")){
			kernel_mat <- as.matrix(kernel_mat)
		}# end fi
		
		if(ncol(kernel_mat) != nrow(kernel_mat)){stop("kernel_mat is a square matrix")}
		num_cell <- nrow(kernel_mat)
		# if(class(norm_counts) != "matrix"){
		if(!is(norm_counts, "matrix")){
			norm_counts <- as.matrix(norm_counts)
		}# end fi
		if(nrow(norm_counts) != num_cell){
			if(ncol(norm_counts) != num_cell){
				stop(" the number of cells in norm_counts and kernel_mat should be matched each other")
			}else{
				norm_counts <- t(norm_counts)
			}# end fi	
		}# end fi
		
		
		## if null add a column vector with 1
		if(!is.matrix(covariates)){
			covariates <- as.matrix(covariates)
		}# end fi
		num_cov <- ncol(covariates)
		## calculate Moore-Penrose pseudoinverse for covariates
		Xdagger <- pinv(covariates)
		
		## calculate SKS - apply S on columns of kernel_mat and then on rows. S is (I - XX+)  
		SKS <- kernel_mat - (covariates%*%Xdagger)%*%kernel_mat
		SKS <- SKS - SKS %*% t(Xdagger) %*% t(covariates)
		
		## eigen
		phis <- Re(eigen(SKS, only.values = TRUE)$values)
		phis[phis<zeros_threshold] <- 0.0
	
		
		
		## calculate k = rank(SKS)
		k <- sum(phis > 0)
		
		## calculate q = dim(ker(SKS) & col(S))
		B <- cbind(kernel_mat, covariates)
		q <- num_cell - sum(svd(B)$d > zeros_threshold)
	 	if(q==0){q <- 1}
		## calculate nominators
		nominators <- apply((SKS%*%norm_counts) * norm_counts, 2, sum)
		
		## calculate dnonimators
		denonimators <- apply((norm_counts - (covariates %*% Xdagger %*% norm_counts)) * norm_counts, 2, sum)
		scores <- nominators / denonimators * (num_cell - num_cov)
		
		## define pvlaues
		pvals <- rep(0,length(scores))
		for(i in seq_len(length(scores))){
			lambda <- c(phis[1:k] - scores[i] / (num_cell - num_cov), ones(1, q) * - scores[i] / (num_cell - num_cov))
			pvals[i] <- davies(0, lambda, acc=1e-7)$Qq
			## to aviod exactaly zeros using liu's method
			if(pvals[i]<0){
				pvals[i] <- liu(0, lambda)
			}# end fi
			if(verbose) {cat(paste("SPARK.SCORE::Davies pvalue 1 = ", pvals[i],"\n"))}
		}#end for
		## return results
		res_test <- data.frame(geneid = colnames(norm_counts), sw=pvals, converged=TRUE)
	}# end fi
  
	#######
	object@res_stest <- res_test
	rm(res_test)
	return(object)
}# end function score test for each kernel




#' Computing gene specific weights for combining multiple p-values
#' 
#' @param object SPARK object
#' @param fit.model The model to be fitted, either "poisson" or "gaussian"
#' @param kernel_mat The kernel matrix
#' @param check_positive Check the kernel matrix is positive or not
#' @param verbose Output fitting information
#' @export
#' 
spark.weights <- function(object, 
						  fit.model,
						  kernel_mat = NULL, 
						  check_positive = TRUE, 
						  verbose = TRUE) {
	## The number of cell
	num_cell <- ncol(object@counts)
	
	## number of genes to testing
	num_gene_test <- length(object@res_vc)
	
	## calculate kernel matrix using location information
	## !!the order of kernel matrices must be consistence with that in spark.test function
	if(is.null(kernel_mat)){ 
		## Euclid distance, and compute the range of parameter l
		ED <- as.matrix(dist(object@location[ ,1:2]))
		lrang <- ComputeGaussianPL(ED, compute_distance=FALSE)[3:7]
		## total 10 kernels, i.e., 5 gaussian and 5 periodic kernels
		kernels <- list(10)
		for(ikernel in c(1:5) ){
			## Gaussian kernel
			if(verbose) {cat(paste0("## testing Gaussian kernel: ",ikernel,"...\n"))}
			kernel_mat <- exp(-ED^2/(2*lrang[ikernel]^2))
			if(check_positive){## kernel matrix should be positive definition matrix
				res <- SysMatEigen(kernel_mat)
				kernels[[2*(ikernel-1)+1]] <- res$kernel_mat
				rm(res)
			}## end fi check
			rm(kernel_mat)
			## Periodic kernel
			if(verbose) {cat(paste0("## testing Periodic kernel: ",ikernel,"...\n"))}
			kernel_mat <- cos(2*pi*ED/lrang[ikernel])
			if(check_positive){## kernel matrix should be positive definition matrix
				res <- SysMatEigen(kernel_mat)
				kernels[[2*ikernel]] <- res$kernel_mat
				rm(res)
			}# end fi check
			rm(kernel_mat)
		}# end for ikernel
	}else{
		kernels <- kernel_mat
	}# end fi
	kernels <- c(kernels, list(diag(num_cell)))
	names(kernels) <- paste0("kernel",1:length(kernels))
	
	#========================================
	if(fit.model == "poisson"){
		resid_working <- NULL
		gene_names <- c()
		for(imodel in 1:num_gene_test){
			model1 <- object@res_vc[[imodel]]
			
			if((class(model1) != "try-error") && (!any(is.na(model1$Y)))){
				resid_working <- cbind(resid_working, model1$Y - model1$X%*%model1$coefficients)
				gene_names <- c(gene_names, names(object@res_vc)[imodel])
			}# end fi
		}# end for
		colnames(resid_working) <- as.character(gene_names)
		rownames(resid_working) <- colnames(object@counts)
	}else if(fit.model == "gaussian"){
		if(is.null(object@scaled_counts)){
			object@scaled_counts <- NormalizeVST(object@counts)
		}## fi
		resid_working <- t(object@scaled_counts)
	}# end fi
	
	
	if(verbose) {cat(paste("## estimating the weights for each gene ...\n"))}
	rest <- ComputeWeightsRcpp(resid_working, length(kernels), kernels)
	
	weights <- abs(rest$weights[,-length(kernels)])/apply(abs(rest$weights[,-length(kernels)]), 1, sum)
	colnames(weights) <- paste0(c("GSP","COS"), rep(1:5,each=2))
	rownames(weights) <- colnames(resid_working)
    object@weights <- weights
	# return the results
    return(object)
}# end function score test

#########################################
#             CODE END                  #
#########################################




