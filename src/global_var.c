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

/* 1) include utils macro */
#include "util.h"
/* 2) include global variables */
/* 3) the current file defs */
#include "global_var.h"
/* 4) other files */
#include "tuple.h"

/********************
 * GLOBALVARIABLES  *
 ********************/



long int SUBMATRIXTABLE[NBMATRICES][4] = {	
    /* Match, Transversion, Transition, other */
    {1, -3, -2, -3},/* Old blast like */
    {2, -3, -2, -3},/* Used somewhere */
    {3, -3, -2, -3},/* Unity matrix   */
    {5, -4, -3, -4},/* EDNAFULL like first variant  */
    {5, -4, -2, -4} /* EDNAFULL like second variant */
};

long int INDELSTABLE[NBMATRICES][2] = {
    {-8,  -2},
    {-12, -4},
    {-16, -4},
    {-16, -4},
    {-16, -4}
};



/* [1]  MODIFIED WHEN STARTING PROGRAM
 *========================================*/


/* 1.1 parametres et flags modifi�s par l'utilisateur */

long int    gp_cost_gap_opened     = -16;
long int    gp_cost_gap_continued  = -4;
long int    gp_costs[4]            = {5,-4,-3,-4}; /* match/mismatch/transition/other*/
long int ** gp_substitution_matrix = NULL;
long int    gp_adhoc_matrix = 0;
long int    gp_matrix = 3;
long int    gp_cost_max_substitution_matrix = 0; /* scoring matrix maximal absolute value */


double      gp_k_blast = -1.0;
double      gp_lambda_blast = -1.0;
long int    gp_xdrop = 25;

double      gp_t = 1e-3; /* t :  max repeats ratio */
char *      gp_motifs[MAX_SEED];
long int    gp_nb_seeds;
long int    gp_seeds_span_max;
long int    gp_seeds_span_min;
long int    gp_seeds_span[MAX_SEED];
long int    gp_seeds_bitweight[MAX_SEED];

long int    gp_mutations_percent = 25;
double      gp_mutations = 0.25;
long int    gp_indels_percent = 8;
double      gp_indels = 0.08;
long int    gp_alpha_percent = 5;
double      gp_alpha = 0.05;

long int    gp_hitcriterion = 2;
long int    gp_reverse = 2;
long int    gp_display = 1;
long int    gp_nbmaxlines = 1000000; /* maximal number of alignements in the output */
long int    gp_distdiag = 16;
long int    gp_lowercase = 0;
double      gp_entropy_min = 2.80;
double      gp_expectation_value = 10.00;

/* post-processing */
long int    gp_win_min = 64;
long int    gp_win_max = 64 * KILO;
double      gp_win_mul = 16.0;

/* sort function index */
long int          gp_sortcriterion            = 70;
SortCrit *        gp_sortcriterion_func       = SortCriterionScore;
long int          gp_sortblockscriterion      = 7;
SortBlocksCrit *  gp_sortblockscriterion_func = SortBlocksCriterionTextNumber;


/* statistical variables */
long int    gp_rho_stat = 0;
long int    gp_delta_stat = 0;
long int    gp_border = 0;

/* files */
long int    gp_nbfiles = 0;
char        gp_files[2][16384];
long int    gp_selection_fasta = 0;

/* data tables */
char *      gp_query       = NULL, *gp_text           = NULL;
char *      gp_query_rev   = NULL;
long int    gp_querysize   = 0,     gp_textsize       = 0;

/* first file chunks */
long int    gp_nbchunks_query  = 0;
char **     gp_chunkname_query = NULL;
long int *  gp_chunksize_query = NULL;
long int *  gp_chunkstrt_query = NULL;

/* second file chunks */
long int    gp_nbchunks_text = 0;
char **     gp_chunkname_text = NULL;
long int *  gp_chunksize_text = NULL;
long int *  gp_chunkstrt_text = NULL;


/* statistical variables */
long int    gp_nb_letters[2][4] =   { { 0, 0, 0, 0
                            },{ 0, 0, 0, 0
			    } }; /* number of 'A','T','G','C' for each file */
double **   gp_freq_letters /* [2][4] */ = NULL ;
                           /* frequency of 'A','T','G','C' for each file */
long int    gp_nb_triplets[2][64] = { { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                            },{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
			    } }; /* number of 'AAA','AAT', ... 'CCC' for each file */
double **   gp_freq_background;        /* 4x4 */
double **   gp_freq_tripletbackground; /* 64x64 */

char lookup[32] =
    { 'A', 'C', 'G', 'T',
      'R', 'Y', 'N', 'U',
      'K', 'W', 'S', 'M',
      'V', 'H', 'D', 'B',
      'a', 'c', 'g', 't',
      'r', 'y', 'n', 'u',
      'k', 'w', 's', 'm',
      'v', 'h', 'd', 'b'
    };


long int backcode[32] = {
  0, 1, 2, 3,
  0, 1, 0, 3,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 1, 2, 3,
  0, 1, 0, 3,
  0, 0, 0, 0,
  0, 0, 0, 0
};

char complement[32] = {
   3, 2, 1, 0,
   5, 4, 6, 0,
   9, 8,11,10,
  15,14,13,12,
  19,18,17,16,
  21,20,22,16,
  25,24,27,26,
  29,28,31,30
};

long int unindexable[32] = {
  0, 0, 0, 0,
  0, 0, 1, 0,
  1, 1, 1, 1,
  1, 1, 1, 1,
  0, 0, 0, 0,
  0, 0, 1, 0,
  1, 1, 1, 1,
  1, 1, 1, 1
};

SortCrit* sortcriteria[NBSORTCRITERIA] = {
  SortCriterionScore,
  SortCriterionEntropy,
  SortCriterionMutual,
  SortCriterionScoreWithEntropy,
  SortCriterionQueryBegin,
  SortCriterionTextBegin,
  SortCriterionPercentIdentityAlign,
  SortCriterionPercentIdentityQuery,
  SortCriterionPercentIdentityText,
};

SortBlocksCrit * sortblockscriteria[NBSORTBLOCKSCRITERIA] =
{
  NULL,
  NULL,
  NULL,
  NULL,
  SortBlocksCriterionQueryNumber,
  SortBlocksCriterionTextNumber,
  SortBlocksCriterionQueryNumber,
  SortBlocksCriterionTextNumber,
  SortBlocksCriterionQueryTextNumber
};


#ifdef TRACE
long int * gp_dots = NULL;
#endif

#ifdef MEM_ALLOCATED
unsigned long int gv_mem_allocated     = 0;
unsigned long int gp_max_mem_allocated = 1024*MEGA;
#endif


/* [2]   MODIFIED WHEN RUNNING ...
 *================================*/

FILE * gv_outstream = NULL;

long int    gv_thread_result = 0;

long int  * gv_delta_shift = NULL;

/* statistics on time spent */
#ifdef STATS
time_t gv_time_spent = 0;
#endif
long int gv_chunk_nb       = 0;
long int gv_chunk_nb_end   = 1;
long int gv_thread_num[MAX_QUERY_CHUNK_THREADS] = {0};
MA *     gv_first_MA       = NULL;
MA *     gv_last_MA        = NULL;
