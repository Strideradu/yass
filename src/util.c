/*
 *  YASS 1.14
 *  Copyright (C) 2004-2011
 *  the YASS team
 *  Laurent Noe, Gregory Kucherov, Mikhail Roytberg, 
 *  Steven Corroy, Antoine De Monte, Christophe Valmir.
 *
 *  laurent.noe|<A>|lifl.fr
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the CeCILL License as published by
 *  the CEA-CNRS-INRIA; either version 2 of the License, or (at your
 *  option) any later version, and the GNU General Public License as
 *  published by the Free Software Foundation; either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This software contains code derived from the GNU libavl library.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "util.h"
#include "global_var.h"

/*
 * renvoie (mant)^(puis)
 *
 */


double dpow(double mant, long int puis)
{
    double c_puis = mant;
    double result = 1;
    long int i;

    for (i = 0; i < (long int) ((sizeof(long int)) << 3); i++) {
	if (puis & 0x00000001) {
	    result *= c_puis;
	}
	c_puis *= c_puis;	/* mant^2  mant^4  mant^8 etc ... */
	puis >>= 1;
    }
    return result;
}


long int ipow(long int mant, long int puis)
{
    long int c_puis = mant;
    long int result = 1;
    long int i;

    for (i = 0; i < (long int) ((sizeof(long int)) << 3); i++) {
	if (puis & 0x00000001) {
	    result *= c_puis;
	}
	c_puis *= c_puis;	/* mant^2  mant^4  mant^8 etc ... */
	puis >>= 1;
    }
    return result;
}



long double ldpow(double mant, long int puis)
{
    long double c_puis = mant;
    long double result = 1;
    long int i;

    for (i = 0; i < (long int) ((sizeof(long int)) << 3); i++) {
	if (puis & 0x00000001) {
	    result *= c_puis;
	}
	c_puis *= c_puis;	/* mant^2  mant^4  mant^8 etc ... */
	puis >>= 1;
    }
    return result;
}


/*
 *
 * Binomial(long int n, long int k)
 *
 */

double C(long int n, long int k)
{
    long int i = 0;
    double num = 1;
    double denum = 1;

    for (i = n - k + 1; i <= n; i++)
	num *= (double) i;
    for (i = 2; i <= k; i++)
	denum *= (double) i;

    return num / denum;
}



/*
 *
 * Taille d'un entier positif affich� en d�cimal
 *
 */

long int size_lint(long int i)
{
    long int result = 0;
    if (i < 0)
	result++;

    if (i == 0)
	result = 1;
    else
	while (i != 0) {
	    i /= 10;
	    result++;
	}
    return result;
}


/*
 * tables de dimension i x j (i lignes -> j colonnes)
 */


long int ** lint_directtable(long int i, long int j)
{
  long int ** itable;
  long int *  jtable;
  long int k;

  itable = (long int **) MALLOC((unsigned) (i * sizeof(long int *)) + (i * j * sizeof(long int)));
  ASSERT(itable,int_directtable);
  jtable = (long int *) (itable + i);

  for(k=0 ; k<i; k++)
    itable[k] = jtable + k * j;

  return itable;
}



void lint_free_directtable(long int ** dtable, long int i, long int j)
{
  FREE(dtable,((unsigned) (i * sizeof(long int *)) + (i * j * sizeof(long int))));
}

double ** dbl_directtable(long int i, long int j)
{

  double ** itable;
  double *  jtable;
  long int k;

  itable = (double **) MALLOC((unsigned) (i * sizeof(double *)) + (i * j * sizeof(double)));
  ASSERT(itable,int_directtable);
  jtable = (double *)  MALLOC(i * j * sizeof(double));
  jtable = (double *) (itable + i);

  for(k=0 ; k<i; k++)
    itable[k] = jtable + k * j;

  return itable;
}


void dbl_free_directtable(double ** dtable, long int i, long int j)
{
  FREE(dtable,((unsigned) (i * sizeof(double *)) + (i * j * sizeof(double))));
}
