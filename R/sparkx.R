#' @title The main function that testing multiple kernels with non-parametric framework
#' @description SPARK-X builds on the covariance kernel test framework, identifying genes with spatial expression pattern in large scale patial transcriptomic studies
#' @param count_in: A n x p gene expression matrix (sparseMatrix)
#' @param locus_in: A n x d location matrix 
#' @param X_in: A n x c covariate matrix
#' @param numCores: A integer specifying multiple threads
#' @param option: A description of kernels to be tested, "single", only tests for projection kernels; "mixture", tests for additional ten kernels.
#' @param verbose: A boolean value whether to print details, for debug purpose
#' @return A list of estimated parameters
#' \item{stats}{A n x k matrix of test statistics for all kernels}
#' \item{res_stest}{A n x k matrix of P values for all kernels}
#' \item{res_mtest}{A n x 2 matrix of combined P values and BY-adjusted P values}
#' @export
sparkx <- function(count_in,locus_in,X_in=NULL,numCores=1,option="mixture",verbose=TRUE){
	if(is(count_in,"matrix")){
		raw_count_mat 	<- as(as.matrix(count_in), "sparseMatrix")
	}else if(is(count_in,"vector")){
		raw_count_mat 	<- as(t(as.matrix(count_in)), "sparseMatrix")
	}else if(is(count_in,"sparseMatrix")){
		raw_count_mat 	<- count_in
	}else{
		stop("counts has to be of following forms: vector,matrix or sparseMatrix")
	}
	rm(count_in)
	

	totalcount 		<- as.vector(sp_sums_Rcpp(raw_count_mat))
	keep_cell_idx 	<- which(totalcount!=0)

	count_mat 		<- raw_count_mat[,keep_cell_idx]
	fil_loc 		<- locus_in[keep_cell_idx,]
	rm(raw_count_mat)

	keep_gene_idx   <- which(as.vector(sp_sums_Rcpp(count_mat,rowSums=TRUE))!=0)
	count_mat 		<- count_mat[keep_gene_idx,]

	numGene 		<- nrow(count_mat)
	numCell			<- ncol(count_mat)

	cat(paste("## ===== SPARK-X INPUT INFORMATION ==== \n"))
	cat(paste("## number of total samples:", numCell,"\n"))
	cat(paste("## number of total genes:", numGene,"\n"))
	if(numCores>1){
		cat(paste("## Running with",numCores,"cores \n"))
	}else{
		cat(paste("## Running with single core, may take some time \n"))
	}

	GeneNames 			<- rownames(count_mat)
	if(sum(is.na(GeneNames))>0){GeneNames[is.na(GeneNames)]<- "NAgene"}
	
	sparkx_list <- list()
	icount 	= 0
	if(option=="mixture"){
		cat(paste0("## Testing With Projection Kernel\n"))
		icount = icount + 1
		final_location 			<- fil_loc
		sparkx_list[[icount]] 	<- sparkx.sk(count_mat,final_location,X_mat=X_in,mc_cores=numCores,verbose= verbose)
		rm(final_location)

		for(iker in 1:5){
			icount = icount + 1
			cat(paste0("## Testing With Gaussian Kernel ",iker,"\n"))
			final_location 			<- apply(fil_loc,2,transloc_func,lker=iker,transfunc="gaussian")
			sparkx_list[[icount]] 	<- sparkx.sk(count_mat,final_location,X_mat=X_in,mc_cores=numCores,verbose= verbose)
			rm(final_location)
		}
		
		for(iker in 1:5){
			icount = icount + 1
			cat(paste0("## Testing With Cosine Kernel ",iker,"\n"))
			final_location 			<- apply(fil_loc,2,transloc_func,lker=iker,transfunc="cosine")
			sparkx_list[[icount]] 	<- sparkx.sk(count_mat,final_location,X_mat=X_in,mc_cores=numCores,verbose= verbose)
			rm(final_location)
		}
		names(sparkx_list) <- c("projection",paste0("gaus",1:5),paste0("cos",1:5))
	}else if(option=="single"){
		cat(paste0("## Testing With Projection Kernel\n"))
		icount = icount + 1
		# scaled_location 	<- apply(fil_loc,2,scale)
		final_location 	<- fil_loc

		sparkx_list[[icount]] 	<- sparkx.sk(count_mat,final_location,X_mat=X_in,mc_cores=numCores,verbose= verbose)
		rm(final_location)
		names(sparkx_list) <- c("projection")
	}else{
		stop("option should be one of following: single or mixture")
	}


	allstat 	<- sapply(sparkx_list,function(x){x$stat})
	allpvals 	<- sapply(sparkx_list,function(x){x$pval})

	rownames(allstat) <- rownames(allpvals) <- GeneNames
	comb_pval 	<- apply(allpvals,1,ACAT)
	pBY 		<- p.adjust(comb_pval,method="BY")

	joint_pval 	<- cbind.data.frame(combinedPval=comb_pval,adjustedPval=pBY)

	final_res 	<- list(stats=allstat,res_stest=allpvals,res_mtest=joint_pval)
	return(final_res)
}



#' @title Testing for single kernel with non-parametric framework
#' @param counts: A n x p gene expression matrix (sparseMatrix)
#' @param infomat: A n x d location matrix 
#' @param X_mat: A n x c covariate matrix
#' @param mc_cores: A integer specifying multiple threads
#' @param verbose: A boolean value whether to print details, for debug purpose
#' @return A data frame with stats and p values for all genes
#' @export
sparkx.sk <- function(counts,infomat,X_mat=NULL,mc_cores=1,verbose=TRUE){

	geneName 			<- rownames(counts)
	if(sum(is.na(geneName))>0){geneName[is.na(geneName)]<- "NAgene"}

	if(is.null(X_mat)){	
		Xinfomat 		<- apply(infomat,2,scale,scale=FALSE)
		loc_inv 		<- solve(t(infomat) %*% infomat)
		kmat_first 		<- Xinfomat %*% loc_inv

		LocDim 			<- ncol(infomat)
		Klam 			<- eigen(crossprod(Xinfomat,kmat_first), only.values=T)$values
		EHL 				<- counts%*%Xinfomat
		numCell 			<- nrow(Xinfomat)

		adjust_nominator <- as.vector(sp_sums_Rcpp(counts^2,TRUE))
		vec_stat 		 <- apply(EHL,1,function(x){x%*%loc_inv%*%as.matrix(x)})*numCell/adjust_nominator

		vec_ybar  		<- as.vector(sp_means_Rcpp(counts,TRUE))
		vec_ylam 		<- unlist(mclapply(1:nrow(counts),function(x){1-numCell*vec_ybar[x]^2/adjust_nominator[x]},mc.cores=mc_cores))
		vec_daviesp 	<- unlist(mclapply(1:nrow(counts),function(x){sparkx_pval(x,vec_ylam,Klam,vec_stat)},mc.cores=mc_cores))
		res_sparkx 		<- as.data.frame(cbind(vec_stat,vec_daviesp))
	}else{
		## check if it is fast or not:YES, much fast
		## otherwise, too much memory required
		if(ncol(counts)<30000){
			counts 				<- as.matrix(counts)
		}

		numCell 			<- nrow(infomat)
		XTX_inv 			<- solve(crossprod(X_mat,X_mat))
		Xadjust_mat 		<- crossprod(infomat,X_mat)%*%crossprod(XTX_inv,t(X_mat))
		Xinfomat 			<- infomat - t(Xadjust_mat)
		info_inv 			<- solve(crossprod(infomat,infomat))
		kmat_first 			<- Xinfomat %*% info_inv
		LocDim 				<- ncol(Xinfomat)

		Klam 			<- eigen(crossprod(Xinfomat,kmat_first), only.values=T)$values

		res_sparkx_list 		<- mclapply(X=1:nrow(counts),FUN=sparkx.sksg,
								expmat= counts,
								xmat = X_mat,
								scaleinfo = Xinfomat,
								numDim = LocDim,
								lambda_K = Klam,
								loc_inv = info_inv,
								mc.cores=mc_cores,
								verbose=verbose)
		res_sparkx 			<- as.data.frame(do.call(rbind,res_sparkx_list))
	}



	colnames(res_sparkx) 	<- c("stat","pval")
	rownames(res_sparkx) 	<- geneName

	return(res_sparkx)
}

#' @title Testing for single gene with each kernel with non-parametric framework
#' @param igene: A gene index
#' @param expmat: A n x p gene expression matrix (sparseMatrix)
#' @param xmat: A n x c covariate matrix
#' @param scaleinfo: A n x d matrix with centered coordinates
#' @param numDim: A integer, the number of dimension of coordinates
#' @param lambda_K: A d-vector of eigen values of kernel being tested
#' @param loc_inv: A d x d matrix, the inverse of the inner product of the location matrix
#' @param verbose: A boolean value whether to print details, for debug purpose
#' @return A vector with a test statistic and p value
#' @export

sparkx.sksg <- function(igene,expmat,xmat,scaleinfo,numDim,lambda_K,loc_inv,verbose=TRUE){
	if(verbose){cat("gene",igene,"\n")}
	single_gene 		<- expmat[igene,]
	numCell 			<- length(single_gene)


	XTX_inv 			<- solve(crossprod(xmat,xmat))
	GTX 				<- crossprod(single_gene,xmat)
	Gadjust_mat 		<- GTX%*%tcrossprod(XTX_inv,GTX)
	# adj_nominator 		<- 1/as.vector(crossprod(single_gene,single_gene)-Gadjust_mat)
	adj_nominator 		<- 1/as.vector(crossprod(single_gene, single_gene))
	lambda_G 			<- as.vector(crossprod(single_gene,single_gene)-Gadjust_mat)*adj_nominator


	YHL <- single_gene%*%scaleinfo
	scoredavies <- adj_nominator*(YHL%*%loc_inv%*%t(YHL))*numCell

	Zsort 				<- sort(lambda_G*lambda_K,decreasing=TRUE)
	results_score 		<- try(davies(scoredavies, Zsort))
	if(class(results_score)!="try-error"){
		pout 				<- results_score$Qq
		if(pout<=0){
			pout 			<- liu(scoredavies, Zsort)
		}
	}else{
		pout <- NA
	}

	return(c(scoredavies,pout))
}

#' @title Calculate SPARK-X P-values
#' @param igene: A gene index
#' @param lambda_G: A p-vector of eigen values for all genes
#' @param lambda_K: A d-vector of eigen values of kernel being tested
#' @param allstat: A p-vector of test statistics for all genes
#' @return A vector of p values
#' @export

sparkx_pval <- function(igene,lambda_G,lambda_K,allstat){
	Zsort 				<- sort(lambda_G[igene]*lambda_K,decreasing=TRUE)
	results_score 		<- try(davies(allstat[igene], Zsort))
	if(class(results_score)!="try-error"){
		pout 				<- results_score$Qq
		if(pout<=0){
			pout 			<- liu(allstat[igene], Zsort)
		}
	}else{
		pout <- NA
	}
	return(pout)
}



#' @title Transforming the coordinate 
#' @param coord A n-vector of coordinate
#' @param lker A index of smoothing or periodic parameter
#' @param transfunc A description of coordinate transform function
#' @export
transloc_func <- function(coord,lker,transfunc="gaussian"){
	l  		<- quantile(coord,probs=seq(0.1,0.9,by=0.2))
	if(transfunc=="gaussian"){
		out <- exp(-coord^2/(2*l[lker]^2))
	}else if(transfunc=="cosine"){
		out <- cos(2*pi*coord/l[lker])
	}
	return(out)
}
