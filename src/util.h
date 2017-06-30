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

#ifndef __UTIL_H_
#define __UTIL_H_
#include <stdio.h>

double dpow(double mant, long int puis);
long int ipow(long int mant, long int puis);
double C(long int n, long int k);
long int size_lint(long int i);
void strnsub(char *buffer, char *motif, long int n, char c);
void unit   (long int in, /* out */ long int * p_out , /* out */ char * p_c);
long int  ** lint_directtable(long int i, long int j);
void         lint_free_directtable(long int ** dtable, long int i, long int j);
double **    dbl_directtable(long int i, long int j);
void         dlb_free_directtable(double ** dtable, long int i, long int j);

/*=======================*
 *      Controls         *
 *=======================*/

/* Prefetching (experimental)
 *
 */

/* #define PREFETCH
 */

/* Enable multi-threading
 *  Pro   : faster
 *  Cons  : THREAD_FORWARD_REVERSE doubles the memory used (but gives strong improvement)
 *          THREAD_ASSEMBLE_ALIGN does not usually improve more than one-half
 *          THREAD_QUERY_CHUNK uses "MAX_QUERY_CHUNK_THREADS" threads when the -S 0 parameter is activated
 */

/* #define THREAD_FORWARD_REVERSE
 */

/* #define THREAD_ASSEMBLE_ALIGN
 */

/* #define THREAD_QUERY_CHUNK
 */

#define MAX_QUERY_CHUNK_THREADS    2


/* Enable modulo memopt
 *  Pro   : less memory used
 *  Cons  : can be sometime faster but more frequently a little bit slower
 *
 */
#define MEMOPT

/* Choose AVL or Red Black trees to regroup alignments
 *  Pro : it seems that AVL are better
 */

#define CHOOSEAVLTREE
/*#define CHOOSERBTREE*/



/* Enable the count of memory allocated that is
 * stopped with the -Alloc size parameter
 *  Pro : more secure
 * Cons : does not work correctly with multithreading
 */

/* #define MEM_ALLOCATED
 */


/* Enable compressed keylist (seed hash table)
 *  Pro  : the search step is faster
 *  Cons : the preprocessing step is slower
 */

#define KEYLISTCOMPRESS


/*
 * Give some stats on the final output
 */

#define STATS

/*
 * Give the plot
 */

#define TRACE

/*
 * memory optimisation on assemble algorithm (slower but runnable on larger
 * sequences ...)
 */

#define DBL_MAX 10e16


/* Enable cache
 *  Pro  : -none-
 *  Cons : allocated more memory
 */

/* #define CACHE
 */



/*
 * first nucleotide position is 0 or 1
 */

/*
 * debug output options
 */

/* #define DEBUG         */
/* #define DEBUG_ASSEMBLE */
/* #define DEBUG_ALIGN    */
/* #define DEBUG_REGROUP */
/* #define DEBUG_DISPLAY*/
/* #define DEBUG_SCORE */
/* #define DEBUG_STATS */
/* #define DEBUG_ALIGNMENT_EXTENSION */
/* #define DEBUG_ALIGNMENT_BORDER */
/* #define DEBUG_ALIGNMENT_STRAIT */

#ifndef PACKAGE_NAME
#define PACKAGE_NAME "yass"
#endif
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "1.12"
#endif


#define ALIGN_EVERY_NB_ITER (4096*4)
#define NBL  80			/* number of colums for alignments and progress bar */
#define TREE_MAX_HEIGHT 256	/* red_black/avl_tree weight */
#define NB_SCORES  3            /* number of scores considered as different */
#define NB_LETTERS 4
#define SENSIBILITY_LAMBDA 1e-9
#define NBLOOPS_K 20		/* used to compute K */
#define DELTA_K 1
#define BCOUNT 1
/* #define BCOUNT 0 */
#define NBSORTCRITERIA 9
#define NBSORTBLOCKSCRITERIA 9
#define NBPOTENTIALCRITERIA 10
#define NBMATRICES 5
#define MAX_SEED 32
#define CHUNKNAMESIZE 1024
#define MAX_WINDOW_REGROUP      ( 128 * 1024)


/*=================*
 * General Macros  *
 *=================*/

#define  RULE "-------------------------------------------------------------------------------\n"
#define DRULE "===============================================================================\n"
#define LOOKUP(i) ((char)(lookup[(int)(i)]))
#define FALSE 0
#define TRUE  1
#define INFINITY_INT  (int)(0x3f000000)
#define ABS(a)   ((a)>(0)?(a):(-(a)))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)>(b)?(b):(a))
#define MIN3(a,b,c) MIN((a),MIN((b),(c)))
#define MAX3(a,b,c) MAX((a),MAX((b),(c)))
#define OUTSTREAM (gv_outstream?gv_outstream:stdout)
#define KILO      (1024)
#define MEGA      (KILO*KILO)


/* A(=0) becomes T(=1) et vise versa .
 * C(=2) becomes G(=3) et vise versa .
 */
#define PURIN(base)      (!((base)& (char)(0x01)))
#define PYRIMIDIN(base)  ( ((base)& (char)(0x01)))
#define COMPLEMENT(base) (complement[(int)(base)])


#define BLACK      0
#define RED        1
#define GREEN      2
#define ORANGE     3
#define BLUE       4
#define PINK       5
#define CYAN       6
#define WHITE      7
/* #define COLOR */
#ifdef COLOR
#define SETCOLOR(std,foreground)  fprintf(std,"\033[%im",30+foreground);
#define RESET(std) fprintf(std,"\033[0m\n");
#else
#define SETCOLOR(std,foreground)
#define RESET(std)
#endif


/*
 * ASSERT (Malloc)
 */

#ifdef MEM_ALLOCATED
#define  ASSERT(prt,errmessage) if(!(prt)){                                                           \
    if (gv_mem_allocated >= gp_max_mem_allocated){                                                    \
       fprintf(stderr,"* Premature end : memory limit of %lu bytes reached : %lu\n",                  \
	       gp_max_mem_allocated, gv_mem_allocated);                                               \
       fprintf(stderr,"                     \"" #prt "\", function \"" #errmessage "()\"\n");         \
       fprintf(stderr,"  please use the \" -Alloc <int> \"parameter to increase this limit \n");      \
       fprintf(stderr,"  or compile with \"./configure --with-low-memory\" \n");                      \
       fprintf(stderr,"  you can also decrease the evalue with the \"-E 1e-6\" parameter\n");         \
    } else                                                                                            \
      fprintf(stderr,"* Allocation Error :\"" #prt "\", function \"" #errmessage "()\"\n");           \
    RESET(stdout)                                                                                     \
    RESET(stderr)                                                                                     \
    exit(0);                                                                                          \
}
#else
#define  ASSERT(prt,errmessage) if(!(prt)){                                              \
  fprintf(stderr,"* Allocation Error :\"" #prt "\", function \"" #errmessage "()\"\n");  \
  RESET(stdout)                                                                          \
  RESET(stderr)                                                                          \
  exit(0);                                                                               \
}
#endif


/*
 * Warning/Errors
 */


#define _WARNING(errmessage)  {                         \
	   fprintf(stderr,"* Warning : " #errmessage);  \
           fprintf(stderr," \n");                       \
}

#define _ERROR(errmessage)  {                        \
	   fprintf(stderr,"* Error : " #errmessage); \
           fprintf(stderr," \n");                    \
           exit(0);                                  \
}


/*
 * Malloc/Free
 */

#ifdef MALLOC_PADDING

#define MALLOC_PREPADDING  4
#define MALLOC_POSTPADDING 4

#ifdef MEM_ALLOCATED

#define MALLOC(size)   (gv_mem_allocated + (unsigned long)(size) < gp_max_mem_allocated?((char *)malloc((size)+ MALLOC_PREPADDING + MALLOC_POSTPADDING) + MALLOCPADDING):(NULL));\
                        gv_mem_allocated +=(unsigned long)(size) + MALLOC_PREPADDING + MALLOC_POSTPADDING;

#define FREE(ptr,size) {free((char *)ptr - MALLOC_PREPADDING);gv_mem_allocated -= size + MALLOC_PREPADDING + MALLOC_POSTPADDING;}

#else

#define MALLOC(size)   ((char *)malloc((size) + MALLOC_PREPADDING + MALLOC_POSTPADDING) + MALLOC_PREPADDING)

#define FREE(ptr,size)  free((char *)ptr - MALLOC_PREPADDING);

#endif

#else

#ifdef MEM_ALLOCATED

#define MALLOC(size)   (gv_mem_allocated + (unsigned long)(size) < gp_max_mem_allocated ? malloc(size):(NULL));\
                        gv_mem_allocated += (unsigned long)(size); /* fprintf(stdout,"[%ld:+%d]\n",gv_mem_allocated,size); */
#define FREE(ptr,size) {free(ptr); gv_mem_allocated -= (unsigned long)(size); /* fprintf(stdout,"[%d:-%d]\n",gv_mem_allocated,size); */ }

#else

#define MALLOC(size)   malloc(size)

#define FREE(ptr,size) free(ptr)

#endif

#endif





#endif
