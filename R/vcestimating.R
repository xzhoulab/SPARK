########################################################################################################################
# Package: SPARK
# Version: 1.0.2
# Date   : 2018-10-23
# Modified: 2019-5-20 08:20:20;2019-8-10 19:24:50; 2020-2-8 00:58:38
# Title  : Count-based spatial model for identifying spatially variable genes
# Authors: S.Q. Sun, J.Q. Zhu, and X. Zhou
# Contacts: shiquans@umich.edu and jiaqiang@umich.edu 
#          University of Michigan, Department of Biostatistics
####################################################################################################

#' Fitting the count-based spatial model to estimate the parameters
#' 
#' @param object SPARK object
#' @param covariates The covariates in experiments, i.e. confounding factors/batch effect
#' @param lib_size The read depth for each cell
#' @param fit.maxiter Iteration
#' @param fit.model The model to be either "poisson" or "gaussian"
#' @param fit.tol Tolerance
#' @param num_core The number of core used when fitting the model
#' @param verbose Output fitting information
#' @export
spark.vc <- function(object, 
					covariates=NULL,
					lib_size=NULL, 
					fit.maxiter=500, 
					fit.tol=1e-5,
					fit.model="poisson",
					num_core=1, 
					verbose=FALSE) {
	## extract the data from the slot of object, createobject() function goes first
	if(length(object@counts) == 0) {
        stop("object@counts has not been set. Run CreateSPARKObject() first and then retry.")
    }# end fi
	
	## load necessary R packages
	#require("CompQuadForm")
	#require("doParallel")
	#require("foreach")
	#require("Matrix")
  
	## set number of core in run
	if(num_core > 1){
		if(num_core>detectCores()){warning("SPARK.VC:: the number of cores you're setting is larger than detected cores!");
		num_core = detectCores()}
	}# end fi
	
	## store the multiple cores
	object@num_core <- num_core
	## covariates, i.e., confounding or batch effects
	if(is.null(covariates)){
		num_cov <- 0
	}else{
		# remove the intercept if added by user, later intercept will add automatically
		if(length(unique(covariates[,1])) == 1){
			covariates <- covariates[, -1]
		}# end fi
		covariates <- as.matrix(covariates)
		num_cov <- ncol(covariates)
	}# end fi
	
	
	## number of genes and cells
	num_gene <- nrow(object@counts)
	num_cell <- ncol(object@counts)
	
	cat(paste("## ===== SPARK INPUT INFORMATION ==== \n"))
	cat(paste("## number of total samples: ", num_cell,"\n"))
	cat(paste("## number of total features: ", num_gene,"\n"))
	cat(paste("## number of adjusted covariates: ", num_cov,"\n"))
	
	object@fit_model <- fit.model
	## main functions
	if(fit.model=="poisson"){
	#*************************************************#
	#   Count-Based Spatial Model Under The Null      #
	#*************************************************#
		cat("# fitting count-based spatial model under the null hypothesis ... \n")
		if(is.null(lib_size)){
			lib_size <- apply(object@counts, 2, sum)
		}else{# already exists
			lib_size <- as.numeric(lib_size)
		}# end fi
		object@lib_size <- lib_size
		#=====================================
		#=====================================
		## do parallel using foreach function
		registerDoParallel(cores = num_core)
		#cl <- makeCluster(num_core)
		#registerDoSNOW(cl)
		#pb <- txtProgressBar(max = num_gene, style = 3)
		#progress <- function(n) setTxtProgressBar(pb, n)
		#opts <- list(progress = progress)
		#res_vc <-foreach(ig = 1:num_gene, .options.snow=opts, .export="spark.fit")%dopar%{
		ig <- 1
		res_vc <- foreach(ig = 1:num_gene)%dopar%{
		#res_vc <- list(num_gene)
		#for(ig in 1:num_gene){
			if(num_cov==0){
				model0 <- try(glm(formula = as.numeric(object@counts[ig,]) ~ 1 + offset(log(lib_size)), family = poisson(link="log")))
				idx <- match(rownames(model.frame(formula = as.numeric(object@counts[ig,]) ~ 1 + offset(log(lib_size)), na.action = na.omit)),rownames(model.frame(formula = as.numeric(object@counts[ig,]) ~ 1 + offset(log(lib_size)), na.action = na.pass)))
			}else{
				model0 <- try(glm(formula = as.numeric(object@counts[ig,]) ~ covariates  + offset(log(lib_size)), family = poisson(link="log")))
				idx <- match(rownames(model.frame(formula = object@counts[ig,] ~ covariates + offset(log(lib_size)), na.action = na.omit)),rownames(model.frame(formula = object@counts[ig,]~covariates + offset(log(lib_size)), na.action = na.pass)))
			}# end fi
	
			## model to fit
			if(verbose) {cat(paste("NO. Gene = ",ig,"\n"))}
			t1 <- system.time(model1 <- try(spark.fit(model0, maxiter = fit.maxiter, tol = fit.tol, verbose=verbose)))

			## store results
			if((class(model1) != "try-error")&&(!any(is.na(model1$Y)))){
				if(verbose){cat(paste("SPARK.CV::tau = ", model1$theta,"\n"))}
				model1$analy_indx <- idx # cell index to run for each gene
				model1$times <- t1[3] # running time for each gene	
			}else{
				model1$converged <- FALSE
			}# end fi
			#######
			model1
			#res_vc[[ig]] <- model1
		########
		}# end for ig, parallel
		#close(pb)
		#stopCluster(cl)
	
		names(res_vc) <- rownames(object@counts)
		object@res_vc <- res_vc
	}else if(fit.model=="gaussian"){
	#************************************************************#
	#    Normalized Count-Based Spatial Model Under The Null     #
	#************************************************************#
		cat("# fitting normalized count-based spatial model under the null hypothesis ... \n")
		## check normalizd counts, vst normalization method
		if(is.null(object@scaled_counts)){
			object@scaled_counts <- NormalizeVST(object@counts)
		}# end fi
		
		## covariates are required in the downstream steps
		## if null add a column vector with 1
		if(is.null(covariates)){
			covariates <- as.matrix(rep(1, num_cell))
		}# end fi
		
		object@res_vc <- list(covariates=as.matrix(covariates))
	}# end fi
	
	# return results
	return(object)
}# end function 



#' Fitting the count-based spatial model to estimate the parameters
#' 
#' @param model0 The initial model fitted by glm function
#' @param maxiter Iteration
#' @param tol Tolerance
#' @param verbose Output fitting information
#' @export
spark.fit <- function(model0, maxiter = 500, tol = 1e-5, verbose = FALSE) {

	if(model0$family$family == "gaussian"){
		tau <- rep(0, 2)
		# fitting vector
		y <- model0$y
		num_cell <- length(y)
		offset <- model0$offset
		if(is.null(offset)) {offset <- rep(0, num_cell)}
	
		#family <- model0$family
		eta <- model0$linear.predictors
		mu  <- model0$fitted.values
		mu.eta <- model0$family$mu.eta(eta)
		D <- mu.eta/sqrt(model0$family$variance(mu))
		tau[2] <- sum(model0$residuals*model0$residuals)/(length(model0$residuals))
		# working vector
		Y <- eta - offset + (y - mu)/mu.eta	#eta is log(y),offset is log(N)
		X <- model.matrix(model0) ## 
		alpha <- model0$coef
		H <- tau[2]*rep(1,num_cell)
		Hinv <- 1.0/H
	
		## number of covariates, including the intercept
		num_cv <- ncol(X)
		
		if(num_cv > 1){## for multiple covariates, time-consuming
			Hinv <- diag(Hinv)
			HinvX <- crossprod(Hinv, X)
			XHinvX <- crossprod(X, HinvX)
			P <- try(Hinv - tcrossprod(tcrossprod(HinvX, chol2inv(chol( XHinvX ))), HinvX))
			Py <- crossprod(P, Y)
			rm(P)
		}else{## only intercept include the model, fit the model more effeciently
			## modified by sun, 2019-8-10 11:04:26
			HinvX <- Hinv
			XHinvX <- sum(HinvX)	
			Py <- Hinv*Y - HinvX%*%(t(HinvX)%*%Y)/XHinvX
		}# end fi num_cv (HinvX%*%t(HinvX)/XHinvX)%*%Y
		model1 <- list(theta = tau, coefficients = alpha, linear.predictors = eta, fitted.values = mu, Y = Y, residuals = y - mu, Py = Py, X = X, D = D)
	}else if(model0$family$family == "poisson"){
		## the number of variance component
		num_vc <- 1	
		fixtau.old <- rep(0, num_vc + 1)
		## to use average information method to fit alternative model
		model1 <- spark.ai(model0, num_vc, maxiter = maxiter, tol = tol, verbose = verbose)
		fixtau.new <- 1*(model1$theta < 1.01 * tol)
		count_num <- 1
		## while(any(fixtau.new != fixtau.old & count_num<50)) {
		while(any(fixtau.new != fixtau.old & count_num<maxiter)) {
			count_num <- count_num + 1
			fixtau.old <- fixtau.new
			## to use average information method to fit alternative model
			model1 <- spark.ai(model0, num_vc, fixtau = fixtau.old, maxiter = maxiter, tol = tol, verbose = verbose)
			fixtau.new <- 1*(model1$theta < 1.01 * tol)
		}# end while
	}# end fi
	
	# return the results
	return(model1)
}# end func


#' Fitting the count-based spatial model using average information algorithm
#' Only for Poisson-based spatial model
#' @param model0 The initial model fitted by glm function
#' @param num_vc The number of random variables
#' @param tau Variance components parameters
#' @param fixtau The parameter to be optimized
#' @param maxiter Iteration
#' @param tol Tolerance
#' @param verbose Output fitting information
#' @export
spark.ai <- function(model0, num_vc, tau = rep(0.1, num_vc+1), fixtau = rep(0, num_vc+1), maxiter = 500, tol = 1e-5, verbose = FALSE){

	# fitting vector
	y <- model0$y

	num_cell <- length(y)
	offset <- model0$offset
	if(is.null(offset)) {offset <- rep(0, num_cell)}
	
	family <- model0$family
	eta <- model0$linear.predictors
	mu  <- model0$fitted.values
	mu.eta <- family$mu.eta(eta)
	D <- mu.eta/sqrt(model0$family$variance(mu))

	# working vector
	Y <- eta - offset + (y - mu)/mu.eta	#eta is log(y),offset is log(N)
	X <- model.matrix(model0) ## 
	alpha <- model0$coef
	
	# number of covariates, including the intercept
	num_cv <- ncol(X)
	# fix first parameter, dispersion
	tau[1] <- 1
	fixtau[1] <- 1

	## find the initial results for tau
	idxtau <- which(fixtau == 0)
	num_cov_mat2 <- sum(fixtau == 0) # is equal to 1
	if(num_cov_mat2 > 0){
		tau[fixtau == 0] <- rep(min(0.9,var(Y)/(num_vc + 1)), num_cov_mat2)
		H <- tau[1]*(1/D^2)
		H <- H + tau[2]
		Hinv <- 1.0/H
	  
		if(num_cv > 1){# for multiple covariates, time-consuming
			Hinv <- diag(Hinv)
			HinvX <- crossprod(Hinv, X)
			XHinvX <- crossprod(X, HinvX)
			P <- try(Hinv - tcrossprod(tcrossprod(HinvX, chol2inv(chol( XHinvX ))), HinvX))
			Py <- crossprod(P, Y)
			tau0 <- tau
			PAPY <- crossprod(P, Py)
			for(ik in 1:num_cov_mat2){
				##modified by sun, 2019-4-13 19:01:22
				tau[idxtau[ik]] <- max(0, tau0[idxtau[ik]] + tau0[idxtau[ik]]^2 * (crossprod(Y, PAPY) - sum(diag(P)))/num_cell)
			}# end for ik loop
			rm(P)
		}else{# only intercept include the model, fit the model more effeciently
			## modified by sun, 2019-4-13 18:57:06
			HinvX <- Hinv
			XHinvX <- sum(HinvX)	
			Py <- Hinv*Y - Hinv*as.numeric(t(HinvX)%*%Y)/XHinvX
			diagP <- Hinv - HinvX*HinvX/XHinvX
  		
			tau0 <- tau
			# modified by sun, 2019-5-20 17:03:52
			#PAPY <- crossprod(P, Py)
			PAPY <- as.numeric( Hinv*Py - Hinv*as.numeric(t(HinvX)%*%Py)/XHinvX)
			for(ik in 1:num_cov_mat2){
				##modified by sun, 2019-4-13 19:01:22
				tau[idxtau[ik]] <- max(0, tau0[idxtau[ik]] + tau0[idxtau[ik]]^2 * (crossprod(Y, PAPY) - sum(diagP))/num_cell)
			}# end for ik loop
			rm(diagP)
			Py <- t(Py)
		}# end fi num_cv
		rm(H)
		rm(Hinv)
		rm(HinvX)
		#rm(Py)
		rm(PAPY)
	}# end if num_cov_mat2
	
	init_tau <- tau
	# with reasonable initial values
	tau[which(tau>10)] <- 0.5
	tau[which(tau<0)] <- 0.5
	#cat(paste0("initial tau=",tau,"\n"))
	for (iter in seq_len(maxiter)){	
		alpha0 <- alpha
		tau0 <- tau
		if(num_cv > 1){# Cppp code to speed up
		  model1 <- CovariatesAI(Y, X, D^2, tau, fixtau, tol)
		}else{# more faster version becasue X is one vector
		  model1 <- noCovariatesAI(Y, X, D^2, tau, fixtau, tol)
		}# end fi
	
		tau <- as.numeric(model1$tau)
		cov <- as.matrix(model1$cov)
		alpha <- as.numeric(model1$alpha)
		eta <- as.numeric(model1$eta) + offset
		
		mu <- family$linkinv(eta)
		mu.eta <- family$mu.eta(eta)
		D <- mu.eta/sqrt(family$variance(mu))	
		Y <- eta - offset + (y - mu)/mu.eta
		## @@@to avoid the strange values@@@, modified by sun, 2019-5-20 18:59:56
		Y[which(abs(Y)>1e3)] <- median(Y)#@@
		## @@@inverse normal to avoid the outlier@@@
		##@@@Y <- qnorm((rank(Y, na.last="keep", ties.method="random")-0.5)/sum(!is.na(Y)));
		# stop rule
		if(2*max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol){	
			break
		}# end fi
		if(max(tau) > tol^(-2)|any(is.infinite(D))|any(is.infinite(mu))|any(is.infinite(eta)) ){
			###== give the initial values obtained by linear regression ==###
			## modified by sun, 2019-5-20 19:34:10
			y <- model0$y
			family <- model0$family
			eta <- model0$linear.predictors
			mu  <- model0$fitted.values
			mu.eta <- family$mu.eta(eta)
			D <- mu.eta/sqrt(model0$family$variance(mu))
			model1$Py <- Py # obtained from initial value
			##working vector
			Y <- eta - offset + (y - mu)/mu.eta	#eta is log(y),offset is log(N)
			X <- model.matrix(model0) ##
			alpha <- model0$coef
			tau <- init_tau
			##
			iter <- maxiter
			break
		}# end fi
	}# end for

	converged <- ifelse(iter < maxiter, TRUE, FALSE)	
	# return results
	return(list(theta = tau, coefficients = alpha, linear.predictors = eta, fitted.values = mu, Y = Y, residuals = y - mu, Py = model1$Py, cov = cov, X = X, D = D, converged = converged))
}# end function


#########################################
#             CODE END                  #
#########################################
