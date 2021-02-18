#include <fstream>
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
// #include <omp.h>


using namespace std;
using namespace arma;
using namespace Rcpp;


//' Do inverse of sysmetric matrix 
//' @param Min A sysmetric matrix
//' 
//' @return A list
//' 
//' @export
// [[Rcpp::export]]
SEXP SysMatEigen(SEXP Min) {
try {
		arma::mat M = as<mat>(Min);
		arma::vec eigval = zeros<vec>( M.n_rows );
		arma::mat eigvec = zeros<mat>( size(M) );
		eig_sym(eigval, eigvec, M, "dc");
		const uvec idx = find(eigval < 1e-8 );
		arma::vec tmp_value = ones<arma::vec>(idx.n_elem);
		eigval.elem( idx ) = tmp_value * 1e-8;
		arma::mat M1 = eigvec.each_row() % eigval.t();
		M = M1 * eigvec.t();
		// return values
		return List::create(Named("eigval") = eigval, Named("eigvec") = eigvec, Named("kernel_mat") = M);
		//return List::create(Named("eigval") = eigval, Named("eigvec") = eigvec);
	}
	catch (std::exception &ex)
	{
		forward_exception_to_r(ex);
	}
	catch (...)
	{
		::Rf_error("C++ exception (unknown reason)...");
	}
	return R_NilValue;
}// end func


//' Do inverse of sparse sysmetric matrix 
//' @param Min A sysmetric matrix
//' @param num_topin The number of top eigen values
//' 
//' @return A list
//' 
//' @export
// [[Rcpp::export]]
SEXP SparseSysMatEigen(SEXP Min, SEXP num_topin) {
try {
		arma::sp_mat M = as<sp_mat>(Min);
		int num_top = Rcpp::as<int>(num_topin);
		arma::vec eigval = zeros<vec>( num_top );
		arma::mat eigvec = zeros<mat>(M.n_rows, num_top );
		eigs_sym(eigval, eigvec, M, num_top);
		const uvec idx = find(eigval < 1e-8 );
		arma::vec tmp_value = ones<arma::vec>(idx.n_elem);
		eigval.elem( idx ) = tmp_value * 1e-8;
		arma::mat M1 = eigvec.each_row() % eigval.t();
		M = M1 * eigvec.t();
		// return values
		return List::create(Named("eigval") = eigval, Named("eigvec") = eigvec, Named("kernel_mat") = M);
		//return List::create(Named("eigval") = eigval, Named("eigvec") = eigvec);
	}
	catch (std::exception &ex)
	{
		forward_exception_to_r(ex);
	}
	catch (...)
	{
		::Rf_error("C++ exception (unknown reason)...");
	}
	return R_NilValue;
}// end func


//' Variance component estimation with covariates using Average Information algorithm
//' @param Yin Working vector
//' @param Xin Covariate matrix
//' @param Din Weight for each gene
//' @param tauin Initial value for variance component
//' @param fixtauin Variance component to be optimized
//' @param tolin Tolerance
//' 
//' @return A list
//' 
//' @export
// [[Rcpp::export]]
SEXP CovariatesAI(SEXP Yin, SEXP Xin, SEXP Din, SEXP tauin, SEXP fixtauin, SEXP tolin) { /*Average Information*/
	try	{
		arma::vec Y = as<vec>(Yin);
		arma::mat X = as<mat>(Xin);
		arma::vec D = as<vec>(Din);
		arma::vec tau = as<vec>(tauin);
		const uvec fixtau = as<uvec>(fixtauin);
		const int num_cov_mat2 = sum(fixtau == 0);
		const double tol = Rcpp::as<double>(tolin);
		uvec ZERO = (tau < tol);

		const int num_cell = X.n_rows;
		//const int num_cvt = X.n_cols; // if number of column X isnot equal to 1
		arma::vec Hinv(num_cell);
		arma::vec one_vec = ones<vec>(num_cell);

		Hinv = tau(0) * (1.0 / (D + 1e-5));
		Hinv += tau(1) * one_vec;
		Hinv = 1.0 / (Hinv + 1e-5);
		arma::vec HinvY = Hinv % Y;
		arma::mat HinvX = X.each_col() % Hinv;
		arma::mat XtHinvX = X.t() * HinvX;
		arma::mat XtHinvX_inv = inv_sympd(XtHinvX);

		arma::mat P = diagmat(Hinv) - HinvX * XtHinvX_inv * HinvX.t();

		arma::vec alpha = XtHinvX_inv * HinvX.t() * Y;
		arma::vec eta = Y - tau(0) * (HinvY - HinvX * alpha) / D;
		arma::vec PY = P * Y;

		if (num_cov_mat2 > 0) {
			const uvec idxtau = find(fixtau == 0);
			arma::mat AImat(num_cov_mat2, num_cov_mat2); //average information matrix
														 //arma::vec PY = P * Y;
			arma::vec score(num_cov_mat2), PAPY;
			for (int i = 0; i < num_cov_mat2; i++) {
				PAPY = P * PY;
				score(i) = dot(Y, PAPY) - sum(P.diag());
				for (int j = 0; j <= i; j++)	{
					AImat(i, j) = dot(PY, PAPY);
					if (j != i)	{
						AImat(j, i) = AImat(i, j);
					} // end fi
				}	 //end for j
			}		  // end for i

			arma::vec Dtau = solve(AImat, score);
			arma::vec tau0 = tau;

			tau.elem(idxtau) = tau0.elem(idxtau) + Dtau;

			tau.elem(find(ZERO % (tau < tol))).zeros();
			double step = 1.0;
			while (any(tau < 0.0)) {
				step *= 0.5;
				tau.elem(idxtau) = tau0.elem(idxtau) + step * Dtau;
				tau.elem(find(ZERO % (tau < tol))).zeros();
			}
			tau.elem(find(tau < tol)).zeros();
		} // end fi
		// boundary tau 0<= tau <=10
		// tau.elem(find(tau >10.0)).ones();
		// return values
		return List::create(Named("tau") = tau, Named("P") = P, Named("cov") = XtHinvX_inv,	Named("alpha") = alpha, Named("Py") = PY, Named("eta") = eta);
	}
	catch (std::exception &ex)
	{
		forward_exception_to_r(ex);
	}
	catch (...)
	{
		::Rf_error("C++ exception (unknown reason)...");
	}
	return R_NilValue;
} // end funcs


//' Variance component estimation without covariates using Average Information algorithm, float format
//' @param Yin Working vector
//' @param Xin Covariate matrix
//' @param Din Weight for each gene
//' @param tauin Initial value for variance component
//' @param fixtauin Variance component to be optimized
//' @param tolin Tolerance
//' 
//' @return A list
//' 
//' 
//' @export
// [[Rcpp::export]]
SEXP noCovariatesAI(SEXP Yin, SEXP Xin, SEXP Din, SEXP tauin, SEXP fixtauin, SEXP tolin) { /*Average Information*/
	try	{
		arma::vec Y = as<vec>(Yin);
		arma::mat X = as<mat>(Xin);
		arma::vec D = as<vec>(Din);
		arma::vec tau = as<vec>(tauin);
		const uvec fixtau = as<uvec>(fixtauin);
		const int num_cov_mat2 = sum(fixtau == 0);
		const double tol = Rcpp::as<double>(tolin);
		uvec ZERO = (tau < tol);

		const int num_cell = X.n_rows;
		//const int num_cvt = X.n_cols; // only suitable for intercept case

		arma::vec Hinv(num_cell);
		arma::vec one_vec = ones<vec>(num_cell);


		Hinv = tau(0) * (1.0 / (D + 1e-5));
		Hinv += tau(1) * one_vec;
		Hinv = 1.0 / (Hinv + 1e-5);

		arma::vec HinvY = Hinv % Y;

		arma::vec HinvX = Hinv;
		double XtHinvX = sum(HinvX);
		double XtHinvX_inv = 1.0 / XtHinvX;
		arma::vec P_diag = Hinv - (HinvX % HinvX) * XtHinvX_inv;
		double alpha = XtHinvX_inv * dot(HinvX, Y);
		arma::vec eta = Y - tau(0) * (HinvY - HinvX * alpha) / D;

		arma::vec PY = HinvY - HinvX * XtHinvX_inv * (HinvX.t() * Y);


		if (num_cov_mat2 > 0) {
			const uvec idxtau = find(fixtau == 0);
			arma::mat AImat(num_cov_mat2, num_cov_mat2); //average information matrix
														 //arma::vec PY = P * Y;
			arma::vec score(num_cov_mat2);
			for (int i = 0; i < num_cov_mat2; i++) {
				
				arma::vec PAPY = Hinv % PY - HinvX * XtHinvX_inv * (HinvX.t() * PY);
	
				score(i) = dot(Y, PAPY) - sum(P_diag);
				for (int j = 0; j <= i; j++)	{
					AImat(i, j) = dot(PY, PAPY);
					if (j != i)	{
						AImat(j, i) = AImat(i, j);
					} // end fi
				}	 //end for j
			}		  // end for i

			arma::vec Dtau = solve(AImat, score);
			arma::vec tau0 = tau;

			tau.elem(idxtau) = tau0.elem(idxtau) + Dtau;
			tau.elem(find(ZERO % (tau < tol))).zeros();
			double step = 1.0;
			while (any(tau < 0.0)) {
				step *= 0.5;
				tau.elem(idxtau) = tau0.elem(idxtau) + step * Dtau;
				tau.elem(find(ZERO % (tau < tol))).zeros();
			} // end while
			tau.elem(find(tau < tol)).zeros();
		} // end fi
		// boundary tau 0<= tau <=10
		//tau.elem(find(tau >10.0)).ones();
		// return values

		return List::create(Named("tau") = tau, Named("Py") = PY, Named("cov") = XtHinvX_inv,	Named("alpha") = alpha, Named("eta") = eta);
	} // end try
	catch (std::exception &ex)
	{
		forward_exception_to_r(ex);
	}
	catch (...)
	{
		::Rf_error("C++ exception (unknown reason)...");
	}
	return R_NilValue;
} // end funcs




//' Compute the testing quantities without covariates, float format
//' @param yin Working vector
//' @param Pyin The vector P*y
//' @param cov_matin Kernel matrix to be tested
//' @param Din Weight for each gene
//' @param tauin Initial value for variance component
//' 
//' @return A list
//' 
//' 
//' @export
// [[Rcpp::export]]
SEXP ComputeTestQuantRcpp_nocov(SEXP yin, SEXP Pyin, SEXP cov_matin, SEXP Din, SEXP tauin)
{
	try
	{
		arma::vec y = as<vec>(yin);
		arma::vec Py = as<vec>(Pyin);
		arma::mat cov_mat = as<mat>(cov_matin);
		arma::vec D = as<vec>(Din);
		arma::vec tau = as<vec>(tauin);

		const int num_cell = y.n_elem;
		arma::vec Hinv(num_cell);
		arma::vec one_vec = ones<vec>(num_cell);

		Hinv = tau(0) * (1.0 / (D + 1e-5));
		Hinv += tau(1) * one_vec;
		Hinv = 1.0 / (Hinv + 1e-5); // Hinv is a diagonal matrix
		arma::vec Hinvy = Hinv % y;

		arma::vec HinvX = Hinv;
		double XtHinvX = sum(HinvX);
	
		arma::mat P = - arma::kron(HinvX, HinvX.t())/XtHinvX;
		P.diag() = P.diag() + Hinv;

		arma::rowvec PKp2 = HinvX.t()*cov_mat;

		arma::mat PK = cov_mat.each_col() % HinvX - arma::kron(HinvX, PKp2)/XtHinvX;

		double trace_PKP = accu(PK % P);
		double newInfoM_p1 = 0.5 * trace(PK * PK);
		double newInfoM = newInfoM_p1 - 0.5 * trace_PKP*trace_PKP/accu(P % P);
		double ee = trace(PK) / 2.0;
		double kk = newInfoM / (2.0 * ee);
		double df = 2.0 * ee * ee / newInfoM;

		arma::vec PKPy = PK * Py;

		double S0 = 0.5 * dot(y, PKPy);
		//cout<<"S0 = " << S0 <<endl;
		// return values
		return List::create(Named("S0") = S0, Named("ee") = ee, Named("infoMp1") = newInfoM_p1, Named("df") = df, Named("kk") = kk);
	} // end try
	catch (std::exception &ex)
	{
		forward_exception_to_r(ex);
	}
	catch (...)
	{
		::Rf_error("C++ exception (unknown reason)...");
	}
	return R_NilValue;
} // end funcs



//' Compute the testing quantities with covariates, float format
//' @param yin Working vector
//' @param Pyin The vector P*y
//' @param Xin Covariate matrix, including the intercept
//' @param cov_matin Kernel matrix to be tested
//' @param Din Weight for each gene
//' @param tauin Initial value for variance component
//' 
//' @return A list
//' 
//' 
//' @export
// [[Rcpp::export]]
SEXP ComputeTestQuantRcpp_cov(SEXP yin, SEXP Pyin, SEXP Xin, SEXP cov_matin, SEXP Din, SEXP tauin)
{
  try
  {
    arma::vec y = as<vec>(yin);
    arma::vec Py = as<vec>(Pyin);
    arma::mat cov_mat = as<mat>(cov_matin);
    arma::mat X = as<mat>(Xin);
    arma::vec D = as<vec>(Din);
    arma::vec tau = as<vec>(tauin);
    
    const int num_cell = y.n_elem;
    arma::vec Hinv(num_cell);
    arma::vec one_vec = ones<vec>(num_cell);
    
    Hinv = tau(0) * (1.0 / (D + 1e-5));
    Hinv += tau(1) * one_vec;
    Hinv = 1.0 / (Hinv + 1e-5); // Hinv is a diagonal matrix
    arma::vec Hinvy = Hinv % y;
    arma::mat HinvX = X.each_col() % Hinv;
    arma::mat XtHinvX = X.t() * HinvX;
    arma::mat XtHinvX_inv = inv_sympd(XtHinvX);
    arma::mat P = diagmat(Hinv) - HinvX * XtHinvX_inv * HinvX.t();

    // modified by sun, 2019-4-13 16:25:06
    arma::mat PK = P*cov_mat;
    double trace_PKP = accu(PK % P);

    // modified by sun, 2019-4-9 12:26:03
    double newInfoM_p1 = 0.5 * trace(PK * PK);
    double newInfoM = newInfoM_p1 - 0.5 * trace_PKP*trace_PKP/accu(P % P);
    double ee = trace(PK) / 2.0;
    double kk = newInfoM / (2.0 * ee);
    double df = 2.0 * ee * ee / newInfoM;
    arma::vec PKPy = PK * Py;

    double S0 = 0.5 * dot(y, PKPy);
    //double ll = 0.0;
    
    // return values
    return List::create(Named("S0") = S0, Named("ee") = ee, Named("infoMp1") = newInfoM_p1, Named("df") = df, Named("kk") = kk);
  } // end try
  catch (std::exception &ex)
  {
    forward_exception_to_r(ex);
  }
  catch (...)
  {
    ::Rf_error("C++ exception (unknown reason)...");
  }
  return R_NilValue;
} // end funcs


//' Compute the testing quantities for linear mixed model, float format
//' @param yin Working vector
//' @param Pyin The vector P*y
//' @param Xin Covariate matrix, including the intercept
//' @param cov_matin Kernel matrix to be tested
//' @param Din Weight for each gene
//' @param tauin Initial value for variance component
//' 
//' @return A list
//' 
//' 
//' @export
// [[Rcpp::export]]
SEXP ComputeTestQuantRcpp_Gaussian(SEXP yin, SEXP Pyin, SEXP Xin, SEXP cov_matin, SEXP Din, SEXP tauin) {
  try {
    arma::vec y = as<vec>(yin);
    arma::vec Py = as<vec>(Pyin);
    arma::mat cov_mat = as<mat>(cov_matin);
    arma::mat X = as<mat>(Xin);
    arma::vec D = as<vec>(Din);
    arma::vec tau = as<vec>(tauin);
    
    const int num_cell = y.n_elem;
    arma::vec Hinv(num_cell);
    arma::vec one_vec = ones<vec>(num_cell);
    
    Hinv = tau(0) * (1.0 / (D + 1e-5));
    Hinv += tau(1) * one_vec;
    Hinv = 1.0 / (Hinv + 1e-5); // Hinv is a diagonal matrix
    arma::vec Hinvy = Hinv % y;
    arma::mat HinvX = X.each_col() % Hinv;
    arma::mat XtHinvX = X.t() * HinvX;
    arma::mat XtHinvX_inv = inv_sympd(XtHinvX);
    arma::mat P = diagmat(Hinv) - HinvX * XtHinvX_inv * HinvX.t();

    // modified by sun, 2019-8-13 11:17:35
	// compute score, i.e., first derivatives
	arma::vec score = zeros<vec>(2);
    arma::mat PK = P*cov_mat; 
	arma::vec PKPy = PK * Py;
	score(0) = 0.5 * (dot(y, PKPy) - trace(PK));
	score(1) = 0.5 * (dot(Py, Py) - trace(P));
   //double score = 0.5 * (dot(y, PKPy) - trace(PK));

    // modified by sun, 2019-8-13 11:17:58
	// compute variance, i.e., second derivatives
	arma::mat Var = zeros<mat>(2, 2);
	Var(0, 0) = 0.5*dot(y, PK*PKPy);
	Var(0, 1) = 0.5*dot(y, PK*P*Py);
	Var(1, 0) = 0.5*dot(y, P*PK*Py);
	Var(1, 1) = 0.5*dot(y, P*P*Py);

	//Var(0, 0) = dot(y, PK*PKPy) - 0.5 * trace(cov_mat*P*cov_mat*P);
	//Var(0, 1) = dot(y, PK*P*Py) - 0.5 * trace(cov_mat*P*P);
	//Var(1, 0) = dot(y, P*PK*Py) - 0.5 * trace(P*cov_mat*P);
	//Var(1, 1) = dot(y, P*P*Py) - 0.5 * trace(P*P);
	//double Var = 0.5 * dot(y, PK*PKPy);
	arma::mat Var_inv = pinv(Var);
	//double Var =  dot(y, PK*PKPy) - 0.5 * trace(cov_mat*P*cov_mat*P);

	// score statistics
    double Ts = dot(score, Var_inv*score);
 
    // return values
    return List::create(Named("Ts") = Ts, Named("score") = score, Named("Var") = Var);
  } // end try
  catch (std::exception &ex)
  {
    forward_exception_to_r(ex);
  }
  catch (...)
  {
    ::Rf_error("C++ exception (unknown reason)...");
  }
  return R_NilValue;
} // end funcs


//' Compute the kernel weights using psedu-counts
//' @param Yin Working vector
//' @param num_kernelin The number of kernel matrices,  10
//' @param Phiin The kernel matrices and identity matrix,  10 + 1 
//' 
//' @return A list
//' 
//' @export
// [[Rcpp::export]]
SEXP ComputeWeightsRcpp(SEXP Yin, SEXP num_kernelin, SEXP Phiin) {
	try	{
		arma::mat Y = as<mat>(Yin);
		const int num_kernel = Rcpp::as<int>(num_kernelin); // number of kernels inputs
		const Rcpp::List Phi(Phiin); // kernel matrix
	
		const int num_gene = Y.n_cols;
		const int num_cell = Y.n_rows; // 
	
		//arma::fcube PHI(num_cell, num_gene, num_kernel);

		arma::mat Q = zeros<mat>(num_kernel, num_kernel);
		for(int i=1; i<=num_kernel; i++) {
		    stringstream kinsi;
			kinsi << "kernel" << i;
			for(int j=i; j<=num_kernel; j++){
				stringstream kinsj;
				kinsj << "kernel" << j;
				Q(i-1, j-1) = accu(symmatl(as<mat>(Phi[kinsi.str()])) % symmatl(as<mat>(Phi[kinsj.str()])) );
				if(i != j){
					Q(j-1, i-1) = Q(i-1, j-1);
				}// end fi
			}// end for
		}// end for
		//cout<<"Q="<<Q<<endl;
		arma::vec y = zeros<vec>(num_cell);
		arma::vec x = zeros<vec>(num_kernel);
		arma::mat weights = zeros<mat>(num_kernel, num_gene);
		//arma::mat Q_inv = inv_sympd( Q );
		arma::mat Q_inv = pinv( Q );
		for(int ig=0; ig<num_gene; ig++) {
			y = Y.col(ig);
			for(int i=1; i<=num_kernel; i++) {
				stringstream kins;
				kins << "kernel" << i;
				x(i-1) =  dot(y, symmatl(as<mat>(Phi[kins.str()])) * y);
			}// end for
			weights.col(ig) = Q_inv * x;
		}// end for

		
		// boundary tau 0<= tau <=10
		// tau.elem(find(tau >10.0)).ones();
		//arma::mat AImat_inv = inv_sympd(AImat);
		// return values
		return List::create(Named("weights") = weights.t());
	}
	catch (std::exception &ex)
	{
		forward_exception_to_r(ex);
	}
	catch (...)
	{
		::Rf_error("C++ exception (unknown reason)...");
	}
	return R_NilValue;
} // end funcs



//' Calculate column-wise or row-wise mean
//' @param sp_data A sparse matrix
//' @param rowMeans A boolean value, whether to calculate row-wise mean
//' @return A n x 1 or p x 1 matrix 
//' @export
// [[Rcpp::export]]
arma::vec sp_means_Rcpp(arma::sp_mat sp_data, bool rowMeans = false) {

  arma::sp_mat norm_col_sums;
  arma::mat tmp_mat;

  if (rowMeans) {
    norm_col_sums = arma::mean(sp_data, 1);
    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_sums.col(0));}
  else {

    norm_col_sums = arma::mean(sp_data, 0);

    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_sums.row(0).t());
  }

  return tmp_mat;
}

//' Calculate column-wise or row-wise sum
//' @param sp_data A sparse matrix
//' @param rowSums A boolean value, whether to calculate row-wise sum
//' @return A n x 1 or p x 1 matrix 
//' @export
// [[Rcpp::export]]
arma::vec sp_sums_Rcpp(arma::sp_mat sp_data, bool rowSums = false) {
  arma::mat tmp_mat;

  arma::sp_mat norm_col_sums;

  if (rowSums) {
    norm_col_sums = arma::sum(sp_data, 1);
    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_sums.col(0));}
  else {
    norm_col_sums = arma::sum(sp_data, 0);
    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_sums.row(0).t());
  }

  return tmp_mat;
}


//' Count column-wise or row-wise non-zero elements of sparse matrix
//' @param sp_data A sparse matrix
//' @param rowNZ A boolean value, whether to count row-wise non-zero elements
//' @return A n x 1 or p x 1 matrix 
//' @export
// [[Rcpp::export]]
arma::vec sp_nz_count_Rcpp(arma::sp_mat sp_data, bool rowNZ = false) {

  arma::sp_mat sp_data_binary = spones(sp_data);
  arma::mat tmp_mat;
  arma::sp_mat norm_col_sums;

  if (rowNZ) {
    norm_col_sums = arma::sum(sp_data_binary, 1);
    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_sums.col(0));}
  else {
    norm_col_sums = arma::sum(sp_data_binary, 0);
    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_sums.row(0).t());
  }

  return tmp_mat;
}


//' Calculate column-wise or row-wise variance
//' @param sp_data A sparse matrix
//' @param rowVars A boolean value, whether to calculate row-wise variance
//' @return A n x 1 or p x 1 matrix 
//' @export
// [[Rcpp::export]]
arma::vec sp_vars_Rcpp(arma::sp_mat sp_data, bool rowVars = false) {
  arma::sp_mat norm_col_Vars;
  arma::mat tmp_mat;
  if (rowVars) {
    norm_col_Vars = arma::var(sp_data,0, 1);
    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_Vars.col(0));}
  else {
    norm_col_Vars = arma::var(sp_data, 0,0);
    tmp_mat = arma::conv_to< arma::mat >::from(norm_col_Vars.row(0).t());
  }


  return tmp_mat;
}

//' Compute Euclidean distance matrix by rows
//' @param x A location matrix
// [[Rcpp::export]]
arma::mat ED_cpp(const arma::mat & x) {
  unsigned int outrows = x.n_rows, i = 0, j = 0;
  double d;
  arma::mat out = zeros<arma::mat>(outrows, outrows);

  for (i = 0; i < outrows - 1; i++) {
    arma::rowvec v1 = x.row(i);
    for (j = i + 1; j < outrows; j++) {
      d = sqrt(sum(pow(v1 - x.row(j), 2.0)));
      out(j, i) = d;
      out(i, j) = d;
    }
  }
  return out;
}


///////////////////////////////////////////////////////////////////////////////////////////
////                             CODE END HERE                                           //
///////////////////////////////////////////////////////////////////////////////////////////
