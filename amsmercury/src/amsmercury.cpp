/*
 * amsmercury.cpp
 * $Id$
 *
 * Copyright (c) 2006 
 * 	Marc Kirchner <marc.kirchner@iwr.uni-heidelberg.de>
 *
 * This code may be freely distributed under the terms of the
 * Lesser GNU Public License (LGPL) version 2 or any later version.
 * 
 */

#include <iostream>
#include <string>
#include <vector>
#include <R.h>
#include <Rdefines.h>
#include "libmercury++.h"

extern "C" {
	SEXP Rmercury(SEXP composition, SEXP charge, SEXP limit);
}

SEXP Rmercury(SEXP composition, SEXP charge, SEXP limit) {

	int ncharge, nlimit, ncomposition;
	int *pcharge, *pcomposition;
	double *plimit;
	
	PROTECT(charge = AS_INTEGER(charge));
	PROTECT(limit = AS_NUMERIC(limit));
	ncharge = LENGTH(charge);
	nlimit = LENGTH(limit);
	if ((ncharge!=1) || (nlimit!=1)) {
		error("mercury: charge and limit must have length==1");
	}
	PROTECT(composition = AS_INTEGER(composition));
	ncomposition = LENGTH(composition);
	pcharge = INTEGER_POINTER(charge);
	plimit = NUMERIC_POINTER(limit);
	pcomposition = INTEGER_POINTER(composition);
	// build std::vector to call mercury::mercury()
	std::vector<unsigned int> comp(ncomposition);
	for (int i = 0; i < ncomposition; i++) {
		comp[i] = pcomposition[i];
	}
	// other: comp(pcomposition[0], pcomposition[0]+ncomposition);
	std::vector<double> mz, ab;
	mercury::mercury(mz, ab, comp, pcharge[0], plimit[0]);

	// build R expression: a list with two vectors
	SEXP result, mzcol, abcol;
	PROTECT(result = NEW_LIST(2));
	PROTECT(mzcol = NEW_NUMERIC(mz.size()));
	double* pmzcol = NUMERIC_POINTER(mzcol);
	// write mz values to vector and save to list
	for (unsigned int j=0; j<mz.size(); j++) {
		pmzcol[j] = mz[j];
	}
	SET_VECTOR_ELT(result, 0, mzcol);
	// write abundance values to vector and save to list
	PROTECT(abcol = NEW_NUMERIC(ab.size()));
	double* pabcol = NUMERIC_POINTER(abcol);
	for (unsigned int j=0; j<ab.size(); j++) {
		pabcol[j] = ab[j];
	}
	SET_VECTOR_ELT(result, 1, abcol);
	// set list names
	SEXP names;
	PROTECT(names = allocVector(STRSXP, 2));
	SET_STRING_ELT(names, 0, mkChar("mz"));
	SET_STRING_ELT(names, 1, mkChar("abundance"));
	SET_NAMES(result, names);
	UNPROTECT(7);
	return result;
}

