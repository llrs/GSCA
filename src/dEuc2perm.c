
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
  Compile: R CMD SHLIB Code/dEuc2perm.c
*/
SEXP dEuc2perm(SEXP cormatrix1, SEXP cormatrix2, SEXP GSdefList, SEXP nP)
{

  /* PREPARATION */
  SEXP GSdef, d;
  int *isAnn, *annIndex, *annIndexPerm, *indexPerm, nAnn = 0, temp = 0;
  int vindex = 0, ng = 0, np = 0, r = 0, c = 0, m = 0, notNA = 0, ri = 0;
  double *sqdvec, *sqdvecPerm, sqd = 0, temp2 = 0;
  time_t rawtime, start, end;

  SEXP cordim = getAttrib(cormatrix1, R_DimSymbol);
  int nGS = LENGTH(GSdefList);
  int G = INTEGER(cordim)[0];
  int P = choose(G, 2);
  int nperm = INTEGER(coerceVector(nP, INTSXP))[0];


  /* CONVERSION */
  // conversion: cormatrix's to a squared difference vector
  PROTECT(cormatrix1 = coerceVector(cormatrix1, REALSXP));
  PROTECT(cormatrix2 = coerceVector(cormatrix2, REALSXP));
  sqdvec = (double *) R_alloc(P, sizeof(double));
  for (int i = 1; i < G; i++){
    for (int j = 0; j < i; j++){

      vindex = G*j - j*(j + 1)/2 + i - j - 1;

      if (!ISNAN(REAL(cormatrix1)[j*G + i]) && 
	  !ISNAN(REAL(cormatrix2)[j*G + i])){
	sqdvec[vindex] = pow(REAL(cormatrix1)[j*G + i] - 
			     REAL(cormatrix2)[j*G + i], 2);
      } else {
	sqdvec[vindex] = NAN;
      };
    };
  };
  UNPROTECT_PTR(cormatrix2);
  UNPROTECT_PTR(cormatrix1);
   /* END OF CONVERSION */



  /* EUCLIDEAN DISTANCE COMPUTATION */
  // distance value to return, of length nGS, along with permutations
  PROTECT(d = allocMatrix(REALSXP, nGS, nperm));

  // to keep track of annotated gene pairs
  isAnn = (int *) R_alloc(P, sizeof(int));
  for (int i = 0; i < P; i++) isAnn[i] = 0; // 0 if not annotated anywhere

  /***** 1. Real data *****/
  // for every gene set
  for (int i = 0; i < nGS; i++){

    // the set-specific correlation vectors 
    PROTECT(GSdef = VECTOR_ELT(GSdefList, i));
    ng = LENGTH(GSdef);
    np = choose(ng, 2);

    // gene-set specific correlation vectors by pair indexing
    sqd = 0; notNA = 0;
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
	temp2 = sqdvec[vindex];
	if (!ISNAN(temp2)){
	  sqd = sqd + temp2; notNA++;
	};
	
	// keep track for permutation
	isAnn[vindex] = 1; // 1 if annotated somewhere

      };
    };
    UNPROTECT_PTR(GSdef);

    if (notNA < np){
      sqd = sqd + (double)(np - notNA) * (sqd/((double)(notNA)));
    };
    REAL(d)[i] = sqrt(sqd)/sqrt((double)np);
  };
  time(&rawtime);
  printf("\n Real data computation completed: %s\n", ctime(&rawtime));




  /***** 2. Permuted data *****/
  // total number of annotated gene pairs
  for (int i = 0; i < P; i++){
    if (isAnn[i] == 1) nAnn++;
  }
  printf(" Starting permuted calculation with\n %i annotated gene pairs out of total %i pairs...\n\n", nAnn, P);

  // index of the annotated gene pairs
  annIndex = (int *) R_alloc(nAnn, sizeof(int));
  m = 0;
  for (int i = 0; i < P; i++){
    if (isAnn[i] == 1){
      annIndex[m] = i; m++;
    };
  };

  // the permutation index vector of annotated gene pairs
  sqdvecPerm = (double *) R_alloc(P, sizeof(double));
  for (int i = 0; i < P; i++) sqdvecPerm[i] = sqdvec[i];
  indexPerm = (int *) R_alloc(P, sizeof(int));
  for (int i = 0; i < P; i++) indexPerm[i] = i;
  annIndexPerm = (int *) R_alloc(nAnn, sizeof(int));
  for (int i = 0; i < nAnn; i++) annIndexPerm[i] = annIndex[i]; // initialize

  // for every permutation
  srand(nperm);
  for (int a = 1; a < nperm; a++){

    start = time(NULL);
    printf(" Permutation %i: ", a + 1);
    // generate the a-th permutation index
    for (int b = 0; b < nAnn; b++){
      ri = (int)((double)nAnn*rand()/((double)RAND_MAX + 1.0));
      temp = annIndexPerm[ri];
      annIndexPerm[ri] = annIndexPerm[b];
      annIndexPerm[b] = temp;
    };

    for (int i = 0; i < nAnn; i++) indexPerm[annIndex[i]] = annIndexPerm[i];
    for (int i = 0; i < nAnn; i++){
      sqdvecPerm[annIndex[i]] = sqdvec[indexPerm[i]];
    };

    // for every gene set
    for (int i = 0; i < nGS; i++){

      // the set-specific correlation vectors 
      PROTECT(GSdef = VECTOR_ELT(GSdefList, i));
      ng = LENGTH(GSdef);
      np = choose(ng, 2);

      // gene-set specific correlation vectors by pair indexing
      sqd = 0; notNA = 0;
      for (int j = 1; j < ng; j++){
	for (int k = 0; k < j; k++){
	
	  // row and column indices of this particular gene pair
	  if (REAL(GSdef)[j] < REAL(GSdef)[k]){
	    r = REAL(GSdef)[k] - 1;
	    c = REAL(GSdef)[j] - 1;
	  } else {
	    r = REAL(GSdef)[j] - 1;
	    c = REAL(GSdef)[k] - 1;
	  };

	  vindex = G*c - c*(c + 1)/2 + r - c - 1;
	  //vindex = indexPerm[vindex];

	  // actual calculation
	  temp2 = sqdvecPerm[vindex];
	  if (!ISNAN(temp2)){
	    sqd = sqd + temp2; notNA++;
	  };

	};
      };
      UNPROTECT_PTR(GSdef);

      if (notNA < np){
	sqd = sqd + (double)(np - notNA) * (sqd/((double)(notNA)));
      };
      REAL(d)[i + a*nGS] = sqrt(sqd)/sqrt((double)np);
    };

    end = time(NULL); time(&rawtime);
    printf("%is, %s", (int)difftime(end, start), ctime(&rawtime));

  };
  /* END OF EUCLIDEAN DISTANCE COMPUTATION */

  UNPROTECT_PTR(d);
  return d;
}


