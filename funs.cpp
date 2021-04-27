//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>
#include <Rmath.h>
#include <omp.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
//[[Rcpp::plugins(openmp)]]

using namespace std;
using namespace Rcpp;
//[[Rcpp::plugins(unwindProtect)]]


// This file habor some functions used in Simulation process 
// 2018-11-22, by julong wei 

//Fun 1, fast linear regression 
// [[Rcpp::export]]
arma::mat lmRcpp(arma::mat & X, arma::mat & Y)
{  
   int nmk = X.n_cols;	
   int nind = X.n_rows;
   int nt = Y.n_cols;
   
   //scale Y
   arma::rowvec ym = arma::mean(Y,0);
   arma::rowvec ysd = arma::stddev(Y,0);
   Y.each_row() -= ym;
   Y.each_row() /= ysd;
   
   //calculate beta
   arma::mat B(nmk, nt);
   for(int i=0; i<nmk; ++i)
   {	
       arma::vec xvec = X.col(i);
       double m = arma::mean(xvec);
       double sd = arma::stddev(xvec);
       arma::vec xs = (xvec-m)/sd;
	   arma::rowvec tmp = (xs.t() * Y)/nind; 
       B.row(i) = tmp;   	
   }
   return B;
}


//Fun 4, expectation 
//[[Rcpp::export()]]
arma::mat calExp(arma::vec & pos, arma::mat & B, arma::mat & C, arma::mat & gen, 
                 const double & LDwindow, const double & cutoff, const int & num_threads )
{
	arma::mat Gmat = gen;
	double w0 = LDwindow;
	double t0 = cutoff;
	int nm = B.n_rows;
	int nindi = Gmat.n_rows;
	
	arma::mat Bexp(nm, 2);
	Bexp.zeros();
	arma::vec Ci = C.col(0) + C.col(1);
	//--------------------------------------//------------------------------------------------------
	//omp_set_num_threads(num_threads);
	//#pragma omp parallel for \
	shared (pos, B, C, Gmat, nm, nindi, w0, t0, Ci, Bexp) \
	default(none) \
	//------------------------------------//------------------------------------------------------
	for (int i=0; i<nm; ++i)
	{		
		double s0 = pos(i) - w0;
		double s1 = pos(i) + w0;		
		arma::mat Csub = C.rows(arma::find((pos>=s0) && (pos<=s1) && (Ci!=0)) );
		if (Csub.n_rows == 0) continue;
		arma::mat Gsub = Gmat.cols(arma::find((pos>=s0) && (pos<=s1) && (Ci!=0)) );
		arma::mat Bsub = B.rows(arma::find((pos>=s0) && (pos<=s1) && (Ci!=0)) );
		arma::mat Bc = Bsub % Csub;
        //scale genotype
		arma::rowvec gm = arma::mean(Gsub,0);
		arma::rowvec gsd = arma::stddev(Gsub,0);
		Gsub.each_row() -= gm;
		Gsub.each_row() /= gsd;
		//		
		arma::vec gi = Gmat.col(i);
		arma::vec Gi = (gi-arma::mean(gi))/arma::stddev(gi);
		arma::vec r = Gsub.t() * Gi * (1.0/nindi);
		arma::vec r2 = pow(r,2);
		arma::uvec sub = arma::find(r2 >= t0);
		arma::rowvec Btmp = arma::trans(r.elem(sub)) * Bc.rows(sub);
		//arma::rowvec Btmp = arma::trans(r) * Bc;
		Bexp.row(i) = Btmp;	
		
	}
	return Bexp;
	
} 


//[[Rcpp::export()]]
arma::mat calLD(arma::vec & pos, arma::mat & B, arma::mat & C, arma::mat & gen, 
                 const double & LDwindow, const double & cutoff, const int & num_threads)
{
	arma::mat Gmat = gen;
	double w0 = LDwindow;
	double t0 = cutoff;
	int nm = B.n_rows;
	int nindi = Gmat.n_rows;
	
	arma::mat Bexp(nm, 2);
	Bexp.zeros();
	arma::vec Ci = C.col(0) + C.col(1);
	
	//--------------------------------------//------------------------------------------------------
//	omp_set_num_threads(num_threads);
//	#pragma omp parallel for \
//	shared (pos, B, C, Gmat, nm, nindi, w0, t0, Ci, Bexp) \
//	default(none) \
	//------------------------------------//------------------------------------------------------
	for (int i=1; i<nm; ++i)
	{		
cout<<i<<endl;
		double s0 = pos(i) - w0;
		double s1 = pos(i) + w0;		
	//	arma::mat Csub = C.rows(arma::find((pos>=s0) && (pos<=s1)));

	//	if (Csub.n_rows == 0) continue;
		arma::mat Gsub = Gmat.cols(arma::find((pos>=s0) && (pos<=s1)));
	//	arma::mat Bsub = B.rows(arma::find((pos>=s0) && (pos<=s1)) );
	//	arma::mat Bc = Bsub % Csub;
        //scale genotype
		
		arma::rowvec gm = arma::mean(Gsub,0);
		arma::rowvec gsd = arma::stddev(Gsub,0);
		Gsub.each_row() -= gm;
		Gsub.each_row() /= gsd;
		//		
		arma::vec gi = Gmat.col(i);
		arma::vec Gi = (gi-arma::mean(gi))/arma::stddev(gi);
		arma::vec r = Gsub.t() * Gi * (1.0/nindi); 
//ofstream ofresult("LD.txt",ios::app);

//ofresult<<setprecision(10)<<r.t();

char file2[10];
sprintf(file2,"%d",i);
char file1[50]="/net/mulan/home/yuef/LD2000/LD_";
char file3[10]=".dat";
strcat(file1,file2);
strcat(file1,file3);
r.save(file1);  
//		arma::vec r2 = pow(r,2);
//		arma::uvec sub = arma::find(r2 >= t0);
//		arma::rowvec Btmp = arma::trans(r.elem(sub)) * Bc.rows(sub);
		//arma::rowvec Btmp = arma::trans(r) * Bc;
//		Bexp.row(i) = Btmp;	
		
	}  
	return Bexp;
	
} 


//[[Rcpp::export()]]
arma::mat calExp2(arma::vec & pos, arma::mat & B, arma::mat & C, 
                 const double & LDwindow, const double & cutoff, const int & num_threads, const int & gene_type)
{
	//arma::mat Gmat = gen;
	double w0 = LDwindow;
	double t0 = cutoff;
	int nm = B.n_rows;
	//int nindi = Gmat.n_rows;
	int type = gene_type;
	arma::mat Bexp(nm, 2);
	Bexp.zeros();
	arma::vec Ci = C.col(0) + C.col(1);
	//--------------------------------------//------------------------------------------------------
	omp_set_num_threads(num_threads);
	#pragma omp parallel for \
	shared (pos, B, C, nm, w0, t0, Ci, Bexp, type) \
	default(none) \
	//------------------------------------//------------------------------------------------------
	for (int i=1; i<nm; ++i)
	{		

		double s0 = pos(i) - w0;
		double s1 = pos(i) + w0;		
		arma::mat Csub = C.rows(arma::find((pos>=s0) && (pos<=s1) && (Ci!=0)) );
		if (Csub.n_rows == 0) continue;
	//	arma::mat Gsub = Gmat.cols(arma::find((pos>=s0) && (pos<=s1) && (Ci!=0)) );
		arma::mat Bsub = B.rows(arma::find((pos>=s0) && (pos<=s1) && (Ci!=0)) );
		arma::mat Bc = Bsub % Csub;
        //scale genotype
	//	arma::rowvec gm = arma::mean(Gsub,0);
	//	arma::rowvec gsd = arma::stddev(Gsub,0);
	//	Gsub.each_row() -= gm;
	//	Gsub.each_row() /= gsd;
		//		
	//	arma::vec gi = Gmat.col(i);
	//	arma::vec Gi = (gi-arma::mean(gi))/arma::stddev(gi);
	//	arma::vec r = Gsub.t() * Gi * (1.0/nindi);
		arma::vec r;
char file2[10];
sprintf(file2,"%d",i);
char file11[50]="/net/mulan/home/yuef/LD2000/LD_";
char file12[50]="/net/mulan/home/yuef/LD2/LD_";
char file13[50]="/net/mulan/home/yuef/LD1/LD_";

char file3[10]=".dat";

if(type==1)
{
strcat(file11,file2);
strcat(file11,file3);
		r.load(file11);
} else if(type==2) {
	strcat(file12,file2);
strcat(file12,file3);
		r.load(file12);
} else {
		strcat(file13,file2);
strcat(file13,file3);
		r.load(file13);
}
        arma::vec Cisub = Csub.col(0) + Csub.col(1); //
		arma::uvec sub2=arma::find(Cisub!=0);//
		r=r.elem(sub2);// for sparse and polygenic setting
		
		arma::vec r2 = pow(r,2);
		arma::uvec sub = arma::find(r2 >= t0);
		arma::rowvec Btmp = arma::trans(r.elem(sub)) * Bc.rows(sub); //out of bounds
	
		//arma::rowvec Btmp = arma::trans(r) * Bc;
		Bexp.row(i) = Btmp;	 
		
	}
	return Bexp;
	
} 

