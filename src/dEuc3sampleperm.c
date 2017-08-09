
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>

/*
  Input: two correlation matrices, a gene set definition list, and
         the desired number of permutations
  Output: a vector of AE distance between two correlation matrices.
  Compile: R CMD SHLIB Code/dEuc2sampleperm.c
*/
SEXP dEuc3sampleperm(SEXP data, SEXP group, SEXP GSdefList, SEXP nP)
{

  /* PREPARATION */
  SEXP GSdef, d;
  int *indexPerm, temp = 0;
  int vindex = 0, ng = 0, np = 0, r = 0, c = 0, ri = 0;
  int notNA12 = 0, notNA13 = 0, notNA23 = 0;
  double *corvec1, *corvec2, *corvec3, *sqdvec12, *sqdvec13, *sqdvec23;
  double *dataPerm, sqd12 = 0, sqd13 = 0, sqd23 = 0, temp2 = 0;
  time_t rawtime, start, end;

  SEXP datadim = getAttrib(data, R_DimSymbol);
  int nGS = LENGTH(GSdefList);
  int G = INTEGER(datadim)[0]; // # genes
  int P = choose(G, 2); // # gene pairs
  int S = INTEGER(datadim)[1]; // # samples
  int n1 = INTEGER(coerceVector(group, INTSXP))[0]; // # samples in group1
  int n2 = INTEGER(coerceVector(group, INTSXP))[1]; // # samples in group2
  int nperm = INTEGER(coerceVector(nP, INTSXP))[0];

  double rho = 0, tempx = 0, tempy = 0;
  double sumx = 0, sumy = 0, sumxy = 0, sumx2 = 0, sumy2 = 0; // for corr calc
  int nn = 0; // # non-NA value in each group (recycled)

  /* PRE-CALCULATION */
  // calculate correlation coefficients and save as a vector for each group
  PROTECT(data = coerceVector(data, REALSXP));
  // group1
  corvec1 = (double *) R_alloc(P, sizeof(double));
  for (int i = 1; i < G; i++){
    for (int j = 0; j < i; j++){
      vindex = G*j - j*(j + 1)/2 + i - j - 1;
      sumx = 0; sumy = 0; sumxy = 0; sumx2 = 0; sumy2 = 0; nn = 0;
      for (int k = 0; k < n1; k++){
	tempx = REAL(data)[k*G + i]; tempy = REAL(data)[k*G + j];
	if (!ISNAN(tempx) && !ISNAN(tempy)){
	  nn++;
	  sumx = sumx + tempx; sumy = sumy + tempy;
	  sumxy = sumxy + tempx*tempy;
	  sumx2 = sumx2 + pow(tempx, 2); sumy2 = sumy2 + pow(tempy, 2);
	};
      };
      rho = ((double)nn * sumxy - sumx * sumy)/
	(sqrt((double)nn * sumx2 - pow(sumx, 2))*
	 sqrt((double)nn * sumy2 - pow(sumy, 2)));
      corvec1[vindex] = rho;
    };
  };
  // group2
  corvec2 = (double *) R_alloc(P, sizeof(double));
  for (int i = 1; i < G; i++){
    for (int j = 0; j < i; j++){
      vindex = G*j - j*(j + 1)/2 + i - j - 1;
      sumx = 0; sumy = 0; sumxy = 0; sumx2 = 0; sumy2 = 0; nn = 0;
      for (int k = n1; k < (n1 + n2); k++){
	tempx = REAL(data)[k*G + i]; tempy = REAL(data)[k*G + j];
	if (!ISNAN(tempx) && !ISNAN(tempy)){
	  nn++;
	  sumx = sumx + tempx; sumy = sumy + tempy;
	  sumxy = sumxy + tempx*tempy;
	  sumx2 = sumx2 + pow(tempx, 2); sumy2 = sumy2 + pow(tempy, 2);
	};
      };
      rho = ((double)nn * sumxy - sumx * sumy)/
	(sqrt((double)nn * sumx2 - pow(sumx, 2))*
	 sqrt((double)nn * sumy2 - pow(sumy, 2)));
      corvec2[vindex] = rho;
    };
  };
  // group3
  corvec3 = (double *) R_alloc(P, sizeof(double));
  for (int i = 1; i < G; i++){
    for (int j = 0; j < i; j++){
      vindex = G*j - j*(j + 1)/2 + i - j - 1;
      sumx = 0; sumy = 0; sumxy = 0; sumx2 = 0; sumy2 = 0; nn = 0;
      for (int k = (n1 + n2); k < S; k++){
	tempx = REAL(data)[k*G + i]; tempy = REAL(data)[k*G + j];
	if (!ISNAN(tempx) && !ISNAN(tempy)){
	  nn++;
	  sumx = sumx + tempx; sumy = sumy + tempy;
	  sumxy = sumxy + tempx*tempy;
	  sumx2 = sumx2 + pow(tempx, 2); sumy2 = sumy2 + pow(tempy, 2);
	};
      };
      rho = ((double)nn * sumxy - sumx * sumy)/
	(sqrt((double)nn * sumx2 - pow(sumx, 2))*
	 sqrt((double)nn * sumy2 - pow(sumy, 2)));
      corvec3[vindex] = rho;
    };
  };

  // calculate pairwise squared differences
  sqdvec12 = (double *) R_alloc(P, sizeof(double));
  sqdvec13 = (double *) R_alloc(P, sizeof(double));
  sqdvec23 = (double *) R_alloc(P, sizeof(double));
  for (int i = 0; i < P; i++){
    sqdvec12[i] = pow(corvec1[i] - corvec2[i], 2);
    sqdvec13[i] = pow(corvec1[i] - corvec3[i], 2);
    sqdvec23[i] = pow(corvec2[i] - corvec3[i], 2);
  };



  /* EUCLIDEAN DISTANCE COMPUTATION */
  // distance value to return, of length nGS, along with permutations
  PROTECT(d = allocMatrix(REALSXP, nGS, nperm));
  /***** 1. Real data *****/
  // for every gene set
  for (int i = 0; i < nGS; i++){

    // set-specific distance
    PROTECT(GSdef = VECTOR_ELT(GSdefList, i));
    ng = LENGTH(GSdef);
    np = choose(ng, 2);
    sqd12 = 0; sqd13 = 0; sqd23 = 0; 
    notNA12 = 0; notNA13 = 0; notNA23 = 0;
    for (int j = 1; j < ng; j++){
      for (int k = 0; k < j; k++){

	// row and column indices of this particular gene pair
	if (REAL(GSdef)[j] < REAL(GSdef)[k]){
	  r = REAL(GSdef)[k] - 1;
	  c = REAL(GSdef)[j] - 1;
	} else {
	  r = REAL(GSdef)[j] - 1;
	  c = REAL(GSdef)[k] - 1;
	}

	vindex = G*c - c*(c + 1)/2 + r - c - 1;
	
	// actual calculation
	// NA handling: impute with the average squared difference
	//              as done by dist() in R
	//              (could ignore the term but will be a smaller value)
	temp2 = sqdvec12[vindex];
	if (!ISNAN(temp2)){
	  sqd12 = sqd12 + temp2; notNA12++;
	};
	temp2 = sqdvec13[vindex];
	if (!ISNAN(temp2)){
	  sqd13 = sqd13 + temp2; notNA13++;
	};
	temp2 = sqdvec23[vindex];
	if (!ISNAN(temp2)){
	  sqd23 = sqd23 + temp2; notNA23++;
	};
	
      };
    };
    UNPROTECT_PTR(GSdef);

    if (notNA12 < np){
      sqd12 = sqd12 + (double)(np - notNA12) * (sqd12/((double)(notNA12)));
    };
    if (notNA13 < np){
      sqd13 = sqd13 + (double)(np - notNA13) * (sqd13/((double)(notNA13)));
    };
    if (notNA23 < np){
      sqd23 = sqd23 + (double)(np - notNA23) * (sqd23/((double)(notNA23)));
    };
    REAL(d)[i] = 
      (sqrt(sqd12) + sqrt(sqd13) + sqrt(sqd23))/3/sqrt((double)np);
  };
  time(&rawtime);
  printf("\n Real data computation completed: %s\n", ctime(&rawtime));



  /***** 2. Permuted data *****/
  // the permuted data matrix
  dataPerm = (double *) R_alloc(G*S, sizeof(double));
  indexPerm = (int *) R_alloc(G*S, sizeof(double));
  for (int i = 0; i < G*S; i++) dataPerm[i] = REAL(data)[i];
  for (int i = 0; i < S; i++) indexPerm[i] = i;

  // set the seed
  srand(nperm);
  
  // for every permutation
  for (int a = 1; a < nperm; a++){

    start = time(NULL);
    printf(" Permutation %i: ", a + 1);

    // generate the a-th permutation index
    for (int b = 0; b < S; b++){
      ri = (int)((double)S*rand()/((double)RAND_MAX + 1.0));
      temp = indexPerm[ri];
      indexPerm[ri] = indexPerm[b];
      indexPerm[b] = temp;
    };

    // assign permuted data into dataPerm
    for (int i = 0; i < S; i++){
      for (int j = 0; j < G; j++){
	dataPerm[i*G + j] = REAL(data)[indexPerm[i]*G + j];
      };
    };

    // permutation-specific correlation vectors
    // group1
    for (int i = 1; i < G; i++){
      for (int j = 0; j < i; j++){
	  vindex = G*j - j*(j + 1)/2 + i - j - 1;
	  sumx = 0; sumy = 0; sumxy = 0; sumx2 = 0; sumy2 = 0; nn = 0;
	  for (int k = 0; k < n1; k++){
	    tempx = dataPerm[k*G + i]; tempy = dataPerm[k*G + j];
	    if (!ISNAN(tempx) && !ISNAN(tempy)){
	      nn++;
	      sumx = sumx + tempx; sumy = sumy + tempy;
	      sumxy = sumxy + tempx*tempy;
	      sumx2 = sumx2 + pow(tempx, 2); sumy2 = sumy2 + pow(tempy, 2);
	    };
	  };
	  rho = ((double)nn * sumxy - sumx * sumy)/
	    (sqrt((double)nn * sumx2 - pow(sumx, 2))*
	     sqrt((double)nn * sumy2 - pow(sumy, 2)));
	  corvec1[vindex] = rho;
      };
    };
    // group2
    for (int i = 1; i < G; i++){
      for (int j = 0; j < i; j++){
	  vindex = G*j - j*(j + 1)/2 + i - j - 1;
	  sumx = 0; sumy = 0; sumxy = 0; sumx2 = 0; sumy2 = 0; nn = 0;
	  for (int k = n1; k < (n1 + n2); k++){
	    tempx = dataPerm[k*G + i]; tempy = dataPerm[k*G + j];
	    if (!ISNAN(tempx) && !ISNAN(tempy)){
	      nn++;
	      sumx = sumx + tempx; sumy = sumy + tempy;
	      sumxy = sumxy + tempx*tempy;
	      sumx2 = sumx2 + pow(tempx, 2); sumy2 = sumy2 + pow(tempy, 2);
	    };
	  };
	  rho = ((double)nn * sumxy - sumx * sumy)/
	    (sqrt((double)nn * sumx2 - pow(sumx, 2))*
	     sqrt((double)nn * sumy2 - pow(sumy, 2)));
	  corvec2[vindex] = rho;
      };
    };
    // group3
    for (int i = 1; i < G; i++){
      for (int j = 0; j < i; j++){
	  vindex = G*j - j*(j + 1)/2 + i - j - 1;
	  sumx = 0; sumy = 0; sumxy = 0; sumx2 = 0; sumy2 = 0; nn = 0;
	  for (int k = (n1 + n2); k < S; k++){
	    tempx = dataPerm[k*G + i]; tempy = dataPerm[k*G + j];
	    if (!ISNAN(tempx) && !ISNAN(tempy)){
	      nn++;
	      sumx = sumx + tempx; sumy = sumy + tempy;
	      sumxy = sumxy + tempx*tempy;
	      sumx2 = sumx2 + pow(tempx, 2); sumy2 = sumy2 + pow(tempy, 2);
	    };
	  };
	  rho = ((double)nn * sumxy - sumx * sumy)/
	    (sqrt((double)nn * sumx2 - pow(sumx, 2))*
	     sqrt((double)nn * sumy2 - pow(sumy, 2)));
	  corvec3[vindex] = rho;
      };
    };

    // calculate pairwise squared differences
    for (int i = 0; i < P; i++){
      sqdvec12[i] = pow(corvec1[i] - corvec2[i], 2);
      sqdvec13[i] = pow(corvec1[i] - corvec3[i], 2);
      sqdvec23[i] = pow(corvec2[i] - corvec3[i], 2);
    };
    


    /* EUCLIDEAN DISTANCE COMPUTATION */
    // for every gene set
    for (int i = 0; i < nGS; i++){
      
      // set-specific distance
      PROTECT(GSdef = VECTOR_ELT(GSdefList, i));
      ng = LENGTH(GSdef);
      np = choose(ng, 2);
      sqd12 = 0; sqd13 = 0; sqd23 = 0; 
      notNA12 = 0; notNA13 = 0; notNA23 = 0;
      for (int j = 1; j < ng; j++){
	for (int k = 0; k < j; k++){
	  
	  // row and column indices of this particular gene pair
	  if (REAL(GSdef)[j] < REAL(GSdef)[k]){
	    r = REAL(GSdef)[k] - 1;
	    c = REAL(GSdef)[j] - 1;
	  } else {
	    r = REAL(GSdef)[j] - 1;
	    c = REAL(GSdef)[k] - 1;
	  }
	  
	  vindex = G*c - c*(c + 1)/2 + r - c - 1;
	  
	  // actual calculation
	  // NA handling: impute with the average squared difference
	  //              as done by dist() in R
	  //              (could ignore the term but will be a smaller value)
	  temp2 = sqdvec12[vindex];
	  if (!ISNAN(temp2)){
	    sqd12 = sqd12 + temp2; notNA12++;
	  };
	  temp2 = sqdvec13[vindex];
	  if (!ISNAN(temp2)){
	    sqd13 = sqd13 + temp2; notNA13++;
	  };
	  temp2 = sqdvec23[vindex];
	  if (!ISNAN(temp2)){
	    sqd23 = sqd23 + temp2; notNA23++;
	  };
	  
	};
      };
      UNPROTECT_PTR(GSdef);
      
      if (notNA12 < np){
	sqd12 = sqd12 + (double)(np - notNA12) * (sqd12/((double)(notNA12)));
      };
      if (notNA13 < np){
	sqd13 = sqd13 + (double)(np - notNA13) * (sqd13/((double)(notNA13)));
      };
      if (notNA23 < np){
	sqd23 = sqd23 + (double)(np - notNA23) * (sqd23/((double)(notNA23)));
      };
      REAL(d)[i + a*nGS] = 
	(sqrt(sqd12) + sqrt(sqd13) + sqrt(sqd23))/3/sqrt((double)np);

    };
    
    end = time(NULL); time(&rawtime);
    printf("%is, %s", (int)difftime(end, start), ctime(&rawtime));
    
  };
  /* END OF EUCLIDEAN DISTANCE COMPUTATION */

  UNPROTECT_PTR(data);
  UNPROTECT_PTR(d);
  return d;
}


