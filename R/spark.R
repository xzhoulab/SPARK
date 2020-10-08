#####################################################################
# Package: SPARK
# Version: 1.0.2
# Date : 2018-10-23
# Modified: 2019-5-20 08:20:20;2019-8-16 09:50:02
# Title : Count-based spatial model for identifying spatially variable genes
# Authors: S.Q. Sun, J.Q. Zhu, and X. Zhou
# Contacts: shiquans@umich.edu and jiaqiang@umich.edu 
#          University of Michigan, Department of Biostatistics
######################################################################
#' Each SPARK object has a number of slots which store information. Key slots to access
#' are listed below.
#'
#' @slot counts The raw expression count matrix
#' @slot scaled_counts Scaled (default is variance stabilizing transformation for each gene) expression matrix; used for visualization
#' @slot weights The weights for combining p-values from multiple kernels
#' @slot kernelmat The kernel matrix for spatial gene recognization, a list
#' @slot lib_size The total read depth for each cell
#' @slot fit_model The model to be fitted, either "poisson" or "gaussian"
#' @slot num_core The number of core used in the package
#' @slot location Cell corrdinates to compute the kernel matrix
#' @slot project Name of the project (for record keeping)
#' @slot res_vc The results for variance component estimation step
#' @slot res_stest The results for variance component testing step for each kernel
#' @slot res_mtest The results for variance component testing step for combined multiple kernels
#' 
setClass("SPARK", slots=list(
  counts = "ANY",
  scaled_counts = "ANY",
  weights = "ANY",
  kernelmat = "list",
  lib_size = "numeric",
  fit_model = "character",
  num_core = "numeric",
  location = "data.frame", 
  project = "character",
  res_vc = "list",
  res_mtest = "data.frame",
  res_stest = "ANY"
) )


#' Create the SPARK object with filtering step
#' @param counts Gene expression count matrix (data.frame), p x n -- p is the number of genes and n is the number of cells
#' @param location Cell location matrix (data.frame) or two-component of t-SNE/UMAP, n x 2
#' @param project Project names
#' @param min_total_counts The minimum counts for each cell for filtering
#' @param percentage The percentage of cells that are expressed for analysis
#' @return Returns SPARK object with filtered gene expression matrix
#' 
#' @export
CreateSPARKObject <- function(counts, location, project = "SPARK", percentage = 0.1, min_total_counts = 10){
  
	## check dimension
	if(ncol(counts)!=nrow(location)){
		stop("The number of cells in counts and location should be consistent! (counts -- p x n; location -- n x 2)")
	}# end fi
  
	## check data order should consistent
	if(!identical(colnames(counts), rownames(location))){
		stop("The column names of counts and row names of location should be should be matched each other! (counts -- p x n; location -- n x 2)")
	}# end fi
  
	### convert into integer
	#counts <- apply(counts, 2, function(x){
	#	x <- ceiling(x)
	#	storage.mode(x) <- 'integer'
	#	return(x) })
  
	## # check data fromat
	#if(!is.data.frame(counts)){
	#  counts <- as.data.frame(counts)
	#}# end fi
  
	## store as sparse matrix
	#if(class(counts) != "dgCMatrix" ){
	#	counts <- as(counts, "dgCMatrix")
	#}# end fi
  
  
	## inheriting
	object <- new(
		Class = "SPARK",
		counts = counts,
		location = location,
		project = project
	)
	rm(counts)
	rm(location)
  
	object.counts <- object@counts
	object.location <- object@location
  
  
	## filtering out lowly expressed genes
	if(percentage > 0){
		gene_use <- which( rowSums(object.counts > 0) >= floor(percentage*ncol(object.counts)) )
		object.counts	<- object.counts[gene_use, ]
	}# end fi
  
	## filter cells that are quite a few number of genes expressed
	cells_use <- colnames(object.counts)
	if(min_total_counts > 0){
		cell_total_counts <- colSums(object.counts)
		cells_use <- cells_use[which(cell_total_counts > min_total_counts)]
		object.counts <- object.counts[, cells_use]
		object.location <- object.location[cells_use, ]
	}# end fi
  
  
	# genes_use <- rownames(object.counts)
	# counts	<- counts[idx,]
	# if(min_cells > 0){
	# 	genes_use <- genes_use[which(rowSums(object.counts > min_counts) > min_cells)]
	# 	object.counts <- object.counts[genes_use, ]
	# }# end fi
  
	object@counts <- object.counts
	object@location <- object.location
	rm(object.counts)
	rm(object.location)
  
	return(object)
}# end function
