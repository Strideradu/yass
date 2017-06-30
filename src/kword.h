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

#ifndef __KWORD_H_
#define __KWORD_H_

#include "threads.h"

/* compute the keyvalue,and the pend, given the pos,data and seed number  */
#define KEY(keyvalue,pend,data,pos,seed) {			            \
    long int  _j_;			                                    \
    pend     = pos;                                                         \
    keyvalue = 0;                                                           \
    for ( _j_ = 0 ; _j_ < gp_seeds_span[seed] ; pend++,_j_++ ){	            \
        long unsigned _data_element_ = (long unsigned)(data[pend]);         \
        char     _seed_element_ = gp_motifs[seed][_j_];                     \
        if  (unindexable[_data_element_]) {                                 \
            keyvalue =  -1;                                                 \
            break;                                                          \
        }                                                                   \
        switch (_seed_element_) {                                           \
           case '#' :                                                       \
              keyvalue <<= 2;                                               \
              keyvalue  |= (backcode[_data_element_]);                      \
           break;                                                           \
           case '@' :                                                       \
              keyvalue <<= 1;                                               \
              keyvalue  |= (backcode[_data_element_]&0x1);                  \
           break;                                                           \
        }                                                                   \
    }                                                                       \
}


/* Read multifasta sequence */
long int CreateData(/* in */  char *filename,
                    /* out */ char **p_data,          /* out */ long int *p_datasize,
                    /* out */ long int *p_nbchunks,   /* out */ char ***p_chunkname,
                    /* out */ long int **p_chunksize, /* out */ long int **p_chunkstart,
                    /* out */ long int * nbletters  /* [4]  */ ,
                    /* out */ long int * nbtriplets /* [64] */ );

/* Sequence reverse complement */
long int CreateReverseComplement(
				 /* in  */ char *   data,         
				 /* in  */ long int datasize, 
				 /* in  */ long int nbchunks,  				     
				 /* in  */ long int * chunksize,
				 /* in  */ long int * chunkstart,
				 /* out */ char ** p_data_rev_out);

/* Compute the number of dots per query chunk */
long int * ComputeDotsTable(long int nbchunks, long int * chunksize);

/* Build the Index */
long int CreateKeyList( /* in */  char *data,         /* in  */ long int datasize,      /* in */ long int seed,
                        /* out */ long int **p_first, /* out */ long int * p_first_size,
                        /* out */ long int **p_list,  /* out */ long int * p_list_size, /* out */ Feature * feature);

long int ComputeLengthAndSortSeeds();

void AutomataOverlaps(long int delph);

long int lowercaseCount(char *data, long int pos, long int len);
long int nCount(char *data, long int pos, long int len);

double Entropy(char *data, long int pos, long int len, long int k);

/* Debugging */
long int DisplaySequence(char *data, long int pos, long int len);
long int DisplayQuery();
long int DisplayText();

long int DisplayData(char * filename, char * data, long int size, long int nbchunks, char ** chunkname, long int * chunksize, long int * chunkstrt);

long int DisplaySeeds();

#endif
