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
#include <string.h>
#include <time.h>

/* 1) include utils macro */
#include "util.h"
/* 2) include global variables */
#include "global_var.h"
/* 3) the current file defs */
#include "assemble.h"

/* 4) other files */
#include "kword.h"
#include "tuple.h"
#include "prdyn.h"
#include "proba.h"
#include "align.h"
#include "display.h"

/*
 * diagonal localy good hash table
 */

#ifdef  CACHE
#define HASH_SHIFT          (7)
#define HASH_DIAG(diagonal) ((diagonal)>>(HASH_SHIFT))
#endif


#define LAST_CHUNK_DIAG_NB  (65535)

/*
 * Modulo access to the last tuple
 */

#ifdef MEMOPT
/* query table index (little more time consuming but memory economic ) */
#define MEMOPT_SETMOD   { long int i=0; memopt_modmask = 0; \
        for (i=0; (                                           \
                   (ALIGN_EVERY_NB_ITER + 1 +                 \
                    keysize_query       + 1 +                 \
                    gp_rho_stat         + 1 +                 \
                    gp_delta_stat       + 1                   \
                    )                                         \
                    >= memopt_modmask)                        \
                  ;i++)                                       \
                 memopt_modmask |=  (1 << i);                 \
  }

#define MEMOPT_MOD      &(memopt_modmask)
#define MEMOPT_TABSIZE   (memopt_modmask+1)
#else
/* full table index (memory consuming) */
#define MEMOPT_SETMOD
#define MEMOPT_MOD
#define MEMOPT_TABSIZE   (keysize_query + keysize_text)
#endif

/*
 * Shift tables used in assemble algorithms
 */


long int initialise_deltashift() {
    long int i, s = 1;
    if (gv_delta_shift == NULL) {
        gv_delta_shift = (long int *) MALLOC ((1 + 2 * gp_delta_stat) * sizeof(long int));
        ASSERT(gv_delta_shift, "initialise_deltashift");

        for (i = 0; i < 1 + 2 * gp_delta_stat; i++) {
            gv_delta_shift[i] = s * (i + 1) / 2;
            s *= -1;
        }
    }
    return 0;
}

/*
 *  definition des fonctions de la table hash ou tableaux
 */
#define INIT_TAB(table, type, size, set)    {                                    \
     table =  (type) MALLOC(size);                                            \
     ASSERT(table,errmessage);                                                \
     memset(table,set,size);                                                  \
}

#define GET_TAB_MIN(table, index, minValue) ((table)[(index) MEMOPT_MOD])
#define GET_TAB(table, index)              ((table)[(index) MEMOPT_MOD])
#define EXIST_TAB(table, index)            ((table)[(index) MEMOPT_MOD])
#define PUT_TAB(table, index, value)        ((table)[(index) MEMOPT_MOD] = value)
#define DEL_TAB(table, index, value)        ((table)[(index) MEMOPT_MOD] = value)
#define FREE_TAB(table, sizeelement)       {FREE(table,(MEMOPT_TABSIZE)*(sizeelement));}
#define RESET_TAB(table, size, set)         {memset(table,set,(size));}


long int Assemble_Single(/*in*/ char *data, /*in*/ long int datasize, Feature *f) {
    fprintf(OUTSTREAM, "Single");
    /* main variables
     * for the linked list computation
     */
    long int *keylist[MAX_SEED];
    long int keylist_size[MAX_SEED];
    long int *first[MAX_SEED];
    long int first_size[MAX_SEED];
    long int keysize_query = datasize;
#ifdef MEMOPT
    long int memopt_modmask = 0;
#else
    long int   keysize_text   = 0;
#endif

    long int *last_tuple_pos_with_diag = NULL;
    tuple **last_tuple_ptr_with_diag = NULL;


#ifdef CACHE
    long int * last_tuple_pos_with_diag_hash = NULL;
#endif


    /* [0] left correction factor set to zero */
    f->left_correction = 0;

    MEMOPT_SETMOD;

    /* to keep lists of tuple found */
    f->last_tl = f->first_tl = CreateTupleList();

    /* [1] key linked list creation */
#ifdef STATS
    f->last_clock = clock();
#endif
    {
        long int seed;
        for (seed = 0; seed < gp_nb_seeds; seed++)
            CreateKeyList(data, datasize, seed, &(first[seed]), &(first_size[seed]), &(keylist[seed]),
                          &(keylist_size[seed]), f);
    }

    /* [2] allocations */
    INIT_TAB(last_tuple_pos_with_diag, long int*, (MEMOPT_TABSIZE) * sizeof(long int), 0xff);
    INIT_TAB(last_tuple_ptr_with_diag, tuple**, (MEMOPT_TABSIZE) * sizeof(tuple *), 0x00);


#ifdef CACHE
    last_tuple_pos_with_diag_hash = (long int *) MALLOC((1 + HASH_DIAG(MEMOPT_TABSIZE)) * sizeof(long int));
    ASSERT(last_tuple_pos_with_diag_hash,"Assemble");
    memset(last_tuple_pos_with_diag_hash,0x7f,(1 + HASH_DIAG(MEMOPT_TABSIZE) * sizeof(long int));
#endif



    /* [3] algorithm to create linked tuples */
    STATS_ADD_CLOCK(f, clock_pre);

    /* for each pos */
    // @nan maybe we just print tuples and should be right for us?
    for (f->i_current = 0; f->i_current < datasize; f->i_current++) {
        /* for each seed */
        long int seed;
        for (seed = 0; seed < gp_nb_seeds; seed++) {

            if (f->i_current <= datasize - gp_seeds_span[seed]) {
                long int index, index_end, i_previous, i_current_end;
                long int code = 0;
                KEY(code, i_current_end, data, f->i_current,
                    seed); /* we take the code of the kword 'w' located at f->i_current position : i_current_end is the first pos "after" the seed */

                if (code >= 0) {

#ifdef KEYLISTCOMPRESS
                    index = first[seed][code];  /* search the first occurrence of 'w' on the sequence */
                    index_end = first[seed][code + 1];
                    i_previous = keylist[seed][index];
#else
                    index      = first[seed][code];  /* search the first occurrence of 'w' on the sequence */
                    index_end  = -1;
                    i_previous = index;
#endif
                    /* consider all the occurrences of 'w' before i_current */
                    while (i_previous < i_current_end - gp_distdiag && index != index_end) {

                        long int diagonal = i_current_end - i_previous;

#ifdef PREFETCH
                        __builtin_prefetch(&GET_TAB(last_tuple_pos_with_diag,diagonal),1,3);
                        __builtin_prefetch(&GET_TAB(last_tuple_ptr_with_diag,diagonal),1,3);
#endif

                        /* consider all the interesting diagonals by adding indel bias  0 +1 -1 +2...
                         * until we reach a new tuple or the upper indel bound */

                        long int observed_diagonal = diagonal;
                        long int delta_diagonal = 0;

#ifdef CACHE
                        {
                          long int i_hash     = HASH_DIAG(MAX((diagonal-gp_delta_stat), 0)MEMOPT_MOD);
                          long int i_hash_end = HASH_DIAG(MIN((diagonal+gp_delta_stat), MEMOPT_TABSIZE)MEMOPT_MOD);
                          do
                            {
                              long int i_previous_hash = last_tuple_pos_with_diag_hash[i_hash];
                              if ((i_previous_hash >= i_current_end - gp_rho_stat) &&
                                  (i_previous_hash <= i_current_end) &&
                                  (i_previous_hash - observed_diagonal <= i_current_end - diagonal) &&
                                  (i_previous_hash >= 0)
                                  )
                                goto diag_start;
                              if (i_hash  == i_hash_end)
                                goto maj_last_tuple_ptr_with_diag;
                              i_hash++;
                              i_hash %= HASH_DIAG(MEMOPT_TABSIZE);
                            }  while (1);

                        }
                      diag_start:
#endif


                        while (delta_diagonal < 2 * gp_delta_stat + 1) {

                            /* search for interesting tuples */
                            long int j = GET_TAB_MIN(last_tuple_pos_with_diag, observed_diagonal,
                                                     i_current_end - gp_rho_stat);

                            /* if tuples found are respecting criteria (neighboors on diagonal and on sequences) */
                            if (j >= 0 &&
                                #ifdef OLDSELECTION
                                j                     <= i_current_end                           &&
                                #endif
                                j >= i_current_end - gp_rho_stat &&
                                #ifdef OLDSELECTION
                                j - observed_diagonal <= i_current_end - diagonal                &&
                                #endif
                                j - observed_diagonal >= i_current_end - diagonal - gp_rho_stat) {

                                tuple *t, *t_;

                                /*---------------------------------------------------------------------*/

                                /* (1) possibly create one instance of tuple [j]
                                 * if this one was not in any list.
                                 */

                                if (EXIST_TAB(last_tuple_ptr_with_diag, observed_diagonal) == 0) {

                                    STATS_NB_CHAINS_BUILT_INC(f);

                                    /* create a tuple [j] with a small seed and update lists */
                                    CREATETUPLE(t, j, observed_diagonal, gp_seeds_span_min);
                                    /* create a tuple list element */
                                    CREATETUPLELIST(f->last_tl->next, t, NULL);
                                    f->last_tl = f->last_tl->next;

                                    PUT_TAB(last_tuple_ptr_with_diag, observed_diagonal, t);
                                } else {

                                    /* (2) we search for the last tuple [j] instanciation
                                     */
                                    t = (tuple *) GET_TAB(last_tuple_ptr_with_diag, observed_diagonal);
                                }

                                /*
                                 * (3) iterate tuples and add this new one to the most adapted one
                                 */

                                do {
                                    long int advance;
                                    if (diagonal == t->diagonal &&
                                        ((advance = i_current_end - t->occurrence) <= gp_seeds_span[seed] + 1)) {

                                        /* (3a) possible fusion with this tuple (no instanciation)
                                         *
                                         * [inter-tuple overlap]
                                         */

                                        if (advance >= 0) {
                                            t->leftsize += advance;
                                            t->occurrence = i_current_end;
                                            goto no_maj_last_tuple_ptr_with_diag;
                                        } else {
                                            goto next_key;
                                        }
                                    }

                                    /* list search for the most fitted */
                                    if (!(t->next))
                                        break;
                                    t = t->next;

                                } while (1);

                                /*
                                 * (3d) any fusion has not been possible
                                 *
                                 * we instanciate the tuple and add this one at the end of the chain.
                                 *
                                 */

                                CREATETUPLE(t_, i_current_end, diagonal, gp_seeds_span[seed]);
                                t->next = t_;
                                PUT_TAB(last_tuple_ptr_with_diag, diagonal, t_);

                                /*
                                 * (4) dont reset last_tuple_ptr_with_diag[diagonal] to NULL
                                 * because it has been modified in (3)
                                 */
                                goto no_maj_last_tuple_ptr_with_diag;

                            }

                            /*---------------------------------------------------------------------*/

                            /* diagonal walk d,d+1,d-1,... */
                            delta_diagonal++;
                            observed_diagonal = diagonal + gv_delta_shift[delta_diagonal];
                            observed_diagonal = MAX(MIN(observed_diagonal, f->i_current), 0);

                            /* border effects avoided*/

                        }/*  while : diagonal walk */


                        if (gp_hitcriterion == 1) {
                            SINGLEHITDIAGONAL(data, datasize, data, datasize);
                        }

                            /* reset last tuple instancied to NULL (because not instance is created)
                             */

#ifdef CACHE
                            maj_last_tuple_ptr_with_diag:
#endif
                        DEL_TAB(last_tuple_ptr_with_diag, diagonal, NULL);

                        no_maj_last_tuple_ptr_with_diag:

                        /*
                         * update last tuple position  with diagonal = 'diagonal'
                         */


                        PUT_TAB(last_tuple_pos_with_diag, diagonal, i_current_end);

#ifdef CACHE
                        last_tuple_pos_with_diag_hash[HASH_DIAG(diagonal MEMOPT_MOD)] = i_current_end;
#endif
                        /* iterate 'w' code list before i_current_end position */
                        next_key:
                        STATS_NB_SEEDS_INC(f);

#ifdef KEYLISTCOMPRESS
                        index++;
                        i_previous = keylist[seed][index];
#else
                        i_previous = index = keylist[seed][index];
#endif
                    } /* while */
                } /* if code >= 0 */
            }
        } /* for each seed*/


        if (!(f->i_current % ALIGN_EVERY_NB_ITER)) {
#ifdef TRACE
            Display_Progress(f->i_current, (keysize_query), f);
#endif

            STATS_ADD_CLOCK(f, clock_chain);

#ifdef DEBUG_ASSEMBLE
            DisplayListTupleList(f->first_tl);
#endif


#ifdef THREAD_ASSEMBLE_ALIGN
            if (f->thread_align) {
              WAIT_THREAD(f->thread_align);
            }
            CREATE_THREAD(f->thread_align,thread_work_align,f);
#else
            AlignAndFree(
                    data, datasize,
                    data, datasize,
                    0,
                    f
            );
#endif

            STATS_ADD_CLOCK(f, clock_align);

        }
    } /* end of [3] */


#ifdef TRACE
    Display_Progress(f->i_current, (keysize_query), f);
#endif

    STATS_ADD_CLOCK(f, clock_chain);

    /* on aligne les tuples restants */
#ifdef DEBUG_ASSEMBLE
    DisplayListTupleList(f->first_tl);
#endif


#ifdef THREAD_ASSEMBLE_ALIGN
    if (f->thread_align) {
      WAIT_THREAD(f->thread_align);
      f->thread_align = 0;
    }
#endif

    /*
     * CREATE_THREAD(f->thread_align,thread_work_align,f);
     */

    AlignAndFree(
            data, datasize,
            data, datasize,
            1,
            f
    );


    STATS_ADD_CLOCK(f, clock_align);

#ifdef CACHE
    FREE(last_tuple_pos_with_diag_hash ,(1 + HASH_DIAG(MEMOPT_TABSIZE)) * sizeof(long int));
#endif

    FREE_TAB(last_tuple_pos_with_diag, sizeof(long int));
    FREE_TAB(last_tuple_ptr_with_diag, sizeof(tuple *));
    {
        long int seed;
        for (seed = 0; seed < gp_nb_seeds; seed++) {
            FREE(keylist[seed], keylist_size[seed]);
            FREE(first[seed], first_size[seed]);
        }
    }
    return 0;
}


long int Assemble_SingleRev(/*in*/   char *datarev, /*in*/  char *data, /*in*/ long int datasize,
                                     Feature *f) {

    fprintf(OUTSTREAM, "Single_Rev");

    /* main variables */
    /* for the linked list computation */
    long int *keylist[MAX_SEED];
    long int keylist_size[MAX_SEED];
    long int *first[MAX_SEED];
    long int first_size[MAX_SEED];
    long int keysize_query = datasize;
#ifndef MEMOPT
    long int   keysize_text   = keysize_query;
#else
    long int memopt_modmask = 0;
#endif

    long int *last_tuple_pos_with_diag = NULL;
    tuple **last_tuple_ptr_with_diag = NULL;

#ifdef CACHE
    long int * last_tuple_pos_with_diag_hash;
#endif


    /* to keep lists of tuple found */
    f->last_tl = f->first_tl = CreateTupleList();

    /* left correction factor set to keysize_query */
    f->left_correction = keysize_query;

    MEMOPT_SETMOD;

    /* [1] key linked list creation */
#ifdef STATS
    f->last_clock = clock();
#endif
    {
        long int seed;
        for (seed = 0; seed < gp_nb_seeds; seed++)
            CreateKeyList(datarev, datasize, seed, &(first[seed]), &(first_size[seed]), &(keylist[seed]),
                          &(keylist_size[seed]), f);
    }

    /* [2] allocations */
    INIT_TAB(last_tuple_pos_with_diag, long int *, (MEMOPT_TABSIZE) * sizeof(long int), 0xff);
    INIT_TAB(last_tuple_ptr_with_diag, tuple **, (MEMOPT_TABSIZE) * sizeof(tuple *), 0x00);

#ifdef CACHE
    last_tuple_pos_with_diag_hash = (long int *) MALLOC((1+HASH_DIAG(MEMOPT_TABSIZE)) * sizeof(long int)));
    ASSERT(last_tuple_pos_with_diag_hash,"Assemble");
    memset(last_tuple_pos_with_diag_hash,0x7f,(1+HASH_DIAG(MEMOPT_TABSIZE) * sizeof(long int)));
#endif



    /* [3] algorithm to create linked tuples */
    STATS_ADD_CLOCK(f, clock_pre);

    /* for each pos */
    for (f->i_current = 0; f->i_current < datasize; f->i_current++) {
        /* for each seed */
        long int seed;
        for (seed = 0; seed < gp_nb_seeds; seed++) {

            if (f->i_current <= datasize - gp_seeds_span[seed]) {
                long int index, index_end, i_previous, i_current_end;
                long int code = 0;
                KEY(code, i_current_end, data, f->i_current,
                    seed); /* we take the code of the kword 'w' located at f->i_current position : i_current_end is the first pos "after" the seed */

                if (code >= 0) {

#ifdef KEYLISTCOMPRESS
                    index = first[seed][code];  /* search the first occurrence of 'w' on the sequence */
                    index_end = first[seed][code + 1];
                    i_previous = keylist[seed][index];
#else
                    index      = first[seed][code];  /* search the first occurrence of 'w' on the sequence */
                    index_end  = -1;
                    i_previous = index;
#endif

                    /* consider all the occurrences of 'w' before f->i_current */


                    while (i_previous < keysize_query - i_current_end && index != index_end) {

                        long int diagonal = i_current_end - i_previous + keysize_query;

#ifdef PREFETCH
                        __builtin_prefetch(&GET_TAB(last_tuple_pos_with_diag,diagonal),1,3);
                        __builtin_prefetch(&GET_TAB(last_tuple_ptr_with_diag,diagonal),1,3);
#endif

                        /* consider all the interesting diagonals by adding indel bias  0 +1 -1 +2...
                         * until we reach a new tuple or the upper indel bound */

                        long int observed_diagonal = diagonal;
                        long int delta_diagonal = 0;

#ifdef CACHE
                        {
                          long int i_hash     = HASH_DIAG(MAX((diagonal-gp_delta_stat), 0)MEMOPT_MOD);
                          long int i_hash_end = HASH_DIAG(MIN((diagonal+gp_delta_stat), MEMOPT_TABSIZE)MEMOPT_MOD);
                          do
                            {
                              long int i_previous_hash = last_tuple_pos_with_diag_hash[i_hash];
                              if ((i_previous_hash >= i_current_end - gp_rho_stat) &&
                                  (i_previous_hash <= i_current_end) &&
                                  (i_previous_hash - observed_diagonal <= i_current_end - diagonal) &&
                                  (i_previous_hash >= 0)
                                  )
                                goto diag_start;
                              if (i_hash  == i_hash_end)
                                goto maj_last_tuple_ptr_with_diag;
                              i_hash++;
                              i_hash %= HASH_DIAG(MEMOPT_TABSIZE);
                            }  while (1);

                        }
                      diag_start:
#endif

                        while (delta_diagonal < 2 * gp_delta_stat + 1) {

                            /* search for interesting tuples */
                            long int j = GET_TAB_MIN(last_tuple_pos_with_diag, observed_diagonal,
                                                     i_current_end - gp_rho_stat);

                            /* if tuples found are respecting criteria (neighboors on diagonal and on sequences) */
                            if (j >= 0 &&
                                #ifdef OLDSELECTION
                                j                     <= i_current_end                           &&
                                #endif
                                j >= i_current_end - gp_rho_stat &&
                                #ifdef OLDSELECTION
                                j - observed_diagonal <= i_current_end - diagonal                &&
                                #endif
                                j - observed_diagonal >= i_current_end - diagonal - gp_rho_stat) {

                                tuple *t, *t_;

                                /*---------------------------------------------------------------------*/

                                /* (1) possibly create one instance of tuple [j]
                                 * if this one was not in any list.
                                 */

                                if (EXIST_TAB(last_tuple_ptr_with_diag, observed_diagonal) == 0) {

                                    STATS_NB_CHAINS_BUILT_INC(f);

                                    /* create a tuple [j] with a small seed and update lists */
                                    CREATETUPLE(t, j, observed_diagonal, gp_seeds_span_min);
                                    /* create a tuple list element */
                                    CREATETUPLELIST(f->last_tl->next, t, NULL);
                                    f->last_tl = f->last_tl->next;

                                    PUT_TAB(last_tuple_ptr_with_diag, observed_diagonal, t);
                                } else {

                                    /* (2) we search for the last tuple [j] instanciation
                                     */
                                    t = (tuple *) GET_TAB(last_tuple_ptr_with_diag, observed_diagonal);
                                }

                                /*
                                 * (3) iterate tuples and add this new one to the most adapted one
                                 */

                                do {
                                    long int advance;
                                    if (diagonal == t->diagonal &&
                                        ((advance = i_current_end - t->occurrence) <= gp_seeds_span[seed] + 1)) {

                                        /* (3a) possible fusion with this tuple (no instanciation)
                                         *
                                         * [inter-tuple overlap]
                                         */

                                        if (advance >= 0) {
                                            t->leftsize += advance;
                                            t->occurrence = i_current_end;
                                            goto no_maj_last_tuple_ptr_with_diag;
                                        } else {
                                            goto next_key;
                                        }
                                    }

                                    /* list search for the most fitted */
                                    if (!(t->next))
                                        break;
                                    t = t->next;

                                } while (1);

                                /*
                                 * (3d) any fusion has not been possible
                                 *
                                 * we instanciate the tuple and add this one at the end of the chain.
                                 *
                                 */

                                CREATETUPLE(t_, i_current_end, diagonal, gp_seeds_span[seed]);
                                t->next = t_;
                                PUT_TAB(last_tuple_ptr_with_diag, diagonal, t_);

                                /*
                                 * (4) dont reset last_tuple_ptr_with_diag[diagonal] to NULL
                                 * because it has been modified in (3)
                                 */
                                goto no_maj_last_tuple_ptr_with_diag;

                            }

                            /*---------------------------------------------------------------------*/

                            /* diagonal walk d,d+1,d-1,... */
                            delta_diagonal++;
                            observed_diagonal = diagonal + gv_delta_shift[delta_diagonal];
                            observed_diagonal = MAX(MIN(observed_diagonal, f->i_current + keysize_query),
                                                    i_current_end);
                            /* border effects avoided*/

                        }/*  while : diagonal walk */


                        if (gp_hitcriterion == 1) {
                            SINGLEHITDIAGONAL(datarev, datasize, data, datasize);
                        }

                            /* reset last tuple instancied to NULL (because not instance is created)
                             */

#ifdef CACHE
                            maj_last_tuple_ptr_with_diag:
#endif
                        DEL_TAB(last_tuple_ptr_with_diag, diagonal, NULL);

                        no_maj_last_tuple_ptr_with_diag:

                        /*
                         * update last tuple position  with diagonal = 'diagonal'
                         */
                        PUT_TAB(last_tuple_pos_with_diag, diagonal, i_current_end);

#ifdef CACHE
                        last_tuple_pos_with_diag_hash[HASH_DIAG(diagonal MEMOPT_MOD)]=i_current_end;
#endif
                        /* iterate 'w' code list before i_current_end position */
                        next_key:
                        STATS_NB_SEEDS_INC(f);

#ifdef KEYLISTCOMPRESS
                        index++;
                        i_previous = keylist[seed][index];
#else
                        i_previous = index = keylist[seed][index];
#endif
                    } /* while i_previous */
                } /* if (code >= 0) */
            }
        } /* for each seed */



        if (!(f->i_current % ALIGN_EVERY_NB_ITER)) {
#ifdef TRACE
            Display_Progress(f->i_current, (keysize_query), f);
#endif

            STATS_ADD_CLOCK(f, clock_chain);

#ifdef DEBUG_ASSEMBLE
            DisplayListTupleList(f->first_tl);
#endif

#ifdef THREAD_ASSEMBLE_ALIGN
            if (f->thread_align) {
              WAIT_THREAD(f->thread_align);
            }
            CREATE_THREAD(f->thread_align,thread_work_align,f);
#else
            AlignAndFree(
                    datarev, datasize,
                    data, datasize,
                    0,
                    f
            );
#endif

            STATS_ADD_CLOCK(f, clock_align);

        }
    }/*end of [3] */


#ifdef TRACE
    Display_Progress(f->i_current, (keysize_query), f);
#endif

    STATS_ADD_CLOCK(f, clock_chain);

    /* on aligne les tuples restants */
#ifdef DEBUG_ASSEMBLE
    DisplayListTupleList(f->first_tl);
#endif


#ifdef THREAD_ASSEMBLE_ALIGN
    if (f->thread_align) {
      WAIT_THREAD(f->thread_align);
      f->thread_align = 0;
    }
#endif

    /*
     * CREATE_THREAD(f->thread_align,thread_work_align,f);
     */

    AlignAndFree(
            datarev, datasize,
            data, datasize,
            1,
            f
    );

    STATS_ADD_CLOCK(f, clock_align);

#ifdef CACHE
    FREE(last_tuple_pos_with_diag_hash,(1 + HASH_DIAG(MEMOPT_TABSIZE) * sizeof(long int));
#endif

    FREE_TAB(last_tuple_pos_with_diag, sizeof(long int));
    FREE_TAB(last_tuple_ptr_with_diag, sizeof(tuple *));
    {
        long int seed;
        for (seed = 0; seed < gp_nb_seeds; seed++) {
            FREE(keylist[seed], keylist_size[seed]);
            FREE(first[seed], first_size[seed]);
        }
    }
    return 0;
}


long int MultiAssemble_Double(char *data_query, long int datasize_query,
                              char *data_text, long int datasize_text,
                              long int nbchunks_text, char **chunkname_text,
                              long int *chunksize_text,
                              long int *chunkstrt_text,
                              Feature *f) {
    //fprintf(OUTSTREAM, "Multiple\n");

    /* used for chunk selection */
    long int i_base = 0;
    long int max_chunk_datasize_text = 0;
    char *chunk_data_text = data_text;

    /*  main variables
     *  for the linked list computation
     */

    long int *keylist_query[MAX_SEED];
    long int keylist_query_size[MAX_SEED];
    long int *first_query[MAX_SEED];
    long int first_query_size[MAX_SEED];

    long int keysize_query = datasize_query;

#ifdef MEMOPT
    long int memopt_modmask = 0;
#else
    long int   keysize_text = 0;
#endif

    long int *last_tuple_pos_with_diag = NULL;
    tuple **last_tuple_ptr_with_diag = NULL;
    /* >> */
    long int *last_chunk_diag_used = NULL;
    long int nb_last_chunk_diag_used = 0;
    /* << */

#ifdef CACHE
    long int *last_tuple_pos_with_diag_hash;
#endif

    /* to keep lists of tuple found */
    f->last_tl = f->first_tl = CreateTupleList();

    {
        long int i_chunk;
        for (i_chunk = 0; i_chunk < nbchunks_text; i_chunk++)
            max_chunk_datasize_text = MAX(max_chunk_datasize_text, chunksize_text[i_chunk]);
    }

#ifndef MEMOPT
    keysize_text = max_chunk_datasize_text;
#endif

    /* [0] left correction factor set to keysize_query */
    f->left_correction = keysize_query;

    MEMOPT_SETMOD;

    /* [1] key linked list creation */
#ifdef STATS
    f->last_clock = clock();
#endif
    {
        long int seed;
        for (seed = 0; seed < gp_nb_seeds; seed++)
            CreateKeyList(data_query, datasize_query, seed, &(first_query[seed]), &(first_query_size[seed]),
                          &(keylist_query[seed]), &(keylist_query_size[seed]), f);
    }

    /* [2] allocations */
    INIT_TAB(last_tuple_pos_with_diag, long int *, (MEMOPT_TABSIZE) * sizeof(long int), 0xff);
    INIT_TAB(last_tuple_ptr_with_diag, tuple **, (MEMOPT_TABSIZE) * sizeof(tuple *), 0x00);

#ifdef CACHE
    last_tuple_pos_with_diag_hash = (long int *) MALLOC((1 + HASH_DIAG(MEMOPT_TABSIZE)) * sizeof(long int));
    ASSERT(last_tuple_pos_with_diag_hash,"Assemble");
    memset(last_tuple_pos_with_diag_hash,0xff,((1 + HASH_DIAG(MEMOPT_TABSIZE)) * sizeof(long int));
#endif

    STATS_ADD_CLOCK(f, clock_pre);

    /* >> */
    last_chunk_diag_used = (long int *) MALLOC(LAST_CHUNK_DIAG_NB * sizeof(long int));
    ASSERT(last_chunk_diag_used, "Assemble");
    memset(last_chunk_diag_used, 0x00, (LAST_CHUNK_DIAG_NB * sizeof(long int)));
    /* << */

    //initialize for iterate tuple
    tuplelist *tl_ = NULL, *tl_prev_ = NULL, *tl_last_;
    tuple *t_ = NULL;

    //volatile long int i_current = f->i_current;
    tl_prev_ = tl_ = f->first_tl;
    tl_last_ = f->last_tl;

    /* for each chunk element */
    for (f->i_chunk = 0; f->i_chunk < nbchunks_text; f->i_chunk++) {

        /*
         * [3] clean used tables  and
         *     and set parameters
         */


        long int chunk_datasize_text = chunksize_text[f->i_chunk];
        long int chunk_keysize_text = chunksize_text[f->i_chunk];
        chunk_data_text = data_text + chunkstrt_text[f->i_chunk];
        i_base = chunkstrt_text[f->i_chunk];

        /* >> */
        if (nb_last_chunk_diag_used < LAST_CHUNK_DIAG_NB) {
            while (nb_last_chunk_diag_used > 0) {
                nb_last_chunk_diag_used--;
                DEL_TAB(last_tuple_pos_with_diag, last_chunk_diag_used[nb_last_chunk_diag_used], -1);
                DEL_TAB(last_tuple_ptr_with_diag, last_chunk_diag_used[nb_last_chunk_diag_used], NULL);
            }
        } else {
            nb_last_chunk_diag_used = 0;
            RESET_TAB(last_tuple_pos_with_diag, (MEMOPT_TABSIZE) * sizeof(long int), 0xff);
            RESET_TAB(last_tuple_ptr_with_diag, (MEMOPT_TABSIZE) * sizeof(tuple *), 0x00);
        }
        /* << */


        /* [4] algorithme de creation des tuples chain�s */
#ifdef DEBUG_ASSEMBLE
        fprintf(stdout,"> f->i_chunk : %ld, chunk_keysize_text = %ld \n",f->i_chunk,chunk_keysize_text);fflush(NULL);
#endif
        /* for each pos  */
        for (f->i_current = 0; f->i_current < chunk_keysize_text; f->i_current++) {
            /* for each seed */
            long int seed;
            for (seed = 0; seed < gp_nb_seeds; seed++) {

                if (f->i_current <= chunk_keysize_text - gp_seeds_span[seed]) {
                    long int index, index_end, i_previous, i_current_end;
                    long int code = 0;
                    KEY(code, i_current_end, chunk_data_text, f->i_current,
                        seed); /* we take the code of the kword 'w' located at f->i_current position : i_current_end is the first pos "after" the seed */

                    if (code >= 0) {

#ifdef KEYLISTCOMPRESS
                        index = first_query[seed][code];  /* search the first occurrence of 'w' on the sequence */
                        index_end = first_query[seed][code + 1];
                        i_previous = keylist_query[seed][index];
#else
                        index      = first_query[seed][code];  /* search the first occurrence of 'w' on the sequence */
                        index_end  = -1;
                        i_previous = index;
#endif

                        /* for each occurrence of 'w' on the "Query" sequence */
                        while (index != index_end) {
                            long int diagonal = i_current_end - i_previous + keysize_query;

#ifdef PREFETCH
                            __builtin_prefetch(&GET_TAB(last_tuple_pos_with_diag,diagonal),1,3);
                            __builtin_prefetch(&GET_TAB(last_tuple_ptr_with_diag,diagonal),1,3);
#endif

                            /* consider all the interesting diagonals by adding indel bias  0 +1 -1 +2...
                             * until we reach a new tuple or the upper indel bound */

                            long int observed_diagonal = diagonal;
                            long int delta_diagonal = 0;

#ifdef CACHE
                            {
                              long int i_hash     = HASH_DIAG(MAX((diagonal-gp_delta_stat), 0)MEMOPT_MOD);
                              long int i_hash_end = HASH_DIAG(MIN((diagonal+gp_delta_stat), MEMOPT_TABSIZE)MEMOPT_MOD);
                              do
                                {
                                  long int i_previous_hash = last_tuple_pos_with_diag_hash[i_hash];
                                  if ((i_previous_hash >= i_current_end - gp_rho_stat) &&
                                      (i_previous_hash <= i_current_end) &&
                                      (i_previous_hash - observed_diagonal <= i_current_end - diagonal) &&
                                      (i_previous_hash >= 0)
                                      )
                                    goto diag_start;
                                  if (i_hash  == i_hash_end)
                                    goto maj_last_tuple_ptr_with_diag;
                                  i_hash++;
                                  i_hash %= HASH_DIAG(MEMOPT_TABSIZE);
                                }  while (1);

                            }
                          diag_start:
#endif

                            while (delta_diagonal < 2 * gp_delta_stat + 1) {

                                /* search for interesting tuples */
                                long int j = GET_TAB_MIN(last_tuple_pos_with_diag, observed_diagonal,
                                                         i_current_end - gp_rho_stat);

                                /* if tuples found are respecting criteria (neighboors on diagonal and on sequences) */
                                if (j >= 0 &&
                                    #ifdef OLDSELECTION
                                    j                     <= i_current_end                           &&
                                    #endif
                                    j >= i_current_end - gp_rho_stat &&
                                    #ifdef OLDSELECTION
                                    j - observed_diagonal <= i_current_end - diagonal                &&
                                    #endif
                                    j - observed_diagonal >= i_current_end - diagonal - gp_rho_stat) {

                                    tuple *t, *t_;

                                    /*---------------------------------------------------------------------*/

                                    /* (1) possibly create one instance of tuple [j]
                                     * if this one was not in any list.
                                     */

                                    if (EXIST_TAB(last_tuple_ptr_with_diag, observed_diagonal) == 0) {

                                        STATS_NB_CHAINS_BUILT_INC(f);

                                        /* create a tuple [j] with a small seed and update lists */
                                        CREATETUPLE(t, j, observed_diagonal, gp_seeds_span_min);
                                        /* create a tuple list element */
                                        CREATETUPLELIST(f->last_tl->next, t, NULL);
                                        f->last_tl = f->last_tl->next;

                                        PUT_TAB(last_tuple_ptr_with_diag, observed_diagonal, t);
                                    } else {

                                        /* (2) we search for the last tuple [j] instanciation
                                         */
                                        t = (tuple *) GET_TAB(last_tuple_ptr_with_diag, observed_diagonal);
                                    }

                                    /*
                                     * (3) iterate tuples and add this new one to the most adapted one
                                     */

                                    do {
                                        long int advance;
                                        if (diagonal == t->diagonal &&
                                            ((advance = i_current_end - t->occurrence) <= gp_seeds_span[seed] + 1)) {

                                            /* (3a) possible fusion with this tuple (no instanciation)
                                             *
                                             * [inter-tuple overlap]
                                             */

                                            if (advance >= 0) {
                                                t->leftsize += advance;
                                                t->occurrence = i_current_end;
                                                goto no_maj_last_tuple_ptr_with_diag;
                                            } else {
                                                goto next_key;
                                            }
                                        }

                                        /* list search for the most fitted */
                                        if (!(t->next))
                                            break;
                                        t = t->next;

                                    } while (1);

                                    /*
                                     * (3d) any fusion has not been possible
                                     *
                                     * we instanciate the tuple and add this one at the end of the chain.
                                     *
                                     */

                                    CREATETUPLE(t_, i_current_end, diagonal, gp_seeds_span[seed]);
                                    t->next = t_;
                                    PUT_TAB(last_tuple_ptr_with_diag, diagonal, t_);

                                    /*
                                     * (4) dont reset last_tuple_ptr_with_diag[diagonal] to NULL
                                     * because it has been modified in (3)
                                     */
                                    goto no_maj_last_tuple_ptr_with_diag;

                                }

                                /*---------------------------------------------------------------------*/

                                /* diagonal walk d,d+1,d-1,... */
                                delta_diagonal++;
                                observed_diagonal = diagonal + gv_delta_shift[delta_diagonal];
                                observed_diagonal = MAX(MIN(observed_diagonal, f->i_current + keysize_query),
                                                        i_current_end);
                                /* border effects avoided*/

                            }/*  while : diagonal walk */


                            if (gp_hitcriterion == 1) {
                                SINGLEHITDIAGONAL_MULTI(data_query, datasize_query, chunk_data_text,
                                                        chunk_datasize_text);
                            }

                                /* reset last tuple instancied to NULL (because not instance is created)
                                 */

#ifdef CACHE
                                maj_last_tuple_ptr_with_diag:
#endif
                            DEL_TAB(last_tuple_ptr_with_diag, diagonal, NULL);

                            no_maj_last_tuple_ptr_with_diag:

                            /*
                             * update last tuple position  with diagonal = 'diagonal'
                             */


                            PUT_TAB(last_tuple_pos_with_diag, diagonal, i_current_end);
                            /* >> */
                            if (nb_last_chunk_diag_used < LAST_CHUNK_DIAG_NB)
                                last_chunk_diag_used[nb_last_chunk_diag_used++] = diagonal;
                            /* << */

#ifdef CACHE
                            last_tuple_pos_with_diag_hash[HASH_DIAG(diagonal MEMOPT_MOD)] = i_current_end;
#endif
                            /* iterate 'w' code list before f->i_current position */
                            next_key:
                            STATS_NB_SEEDS_INC(f);

#ifdef KEYLISTCOMPRESS
                            index++;
                            i_previous = keylist_query[seed][index];
#else
                            i_previous = index = keylist_query[seed][index];
#endif


                        } /* while (index != indexfin) */
                    } /* if (code >= 0) */
                }
            }/* for each seed */



            if (!(f->i_current % ALIGN_EVERY_NB_ITER)) {
#ifdef TRACE
                Display_Progress(f->i_current + i_base, (datasize_text), f);
#endif

                STATS_ADD_CLOCK(f, clock_chain);

#ifdef DEBUG_ASSEMBLE
                DisplayListTupleList(f->first_tl);
#endif


#ifdef THREAD_ASSEMBLE_ALIGN
                if (f->thread_align) {
                  WAIT_THREAD(f->thread_align);
                }
                CREATE_THREAD(f->thread_align,thread_work_align,f);
#else
                AlignAndFree(
                        data_query, datasize_query,
                        chunk_data_text, chunk_datasize_text,
                        0,
                        f);
#endif

                STATS_ADD_CLOCK(f, clock_align);

            }
        }/*end of  [3] (keysize_text loop)  */
#ifdef TRACE
        Display_Progress(f->i_current + i_base, (datasize_text), f);
#endif

        STATS_ADD_CLOCK(f, clock_chain);

#ifdef DEBUG_ASSEMBLE
        DisplayListTupleList(f->first_tl);
#endif


#ifdef THREAD_ASSEMBLE_ALIGN
        if (f->thread_align) {
          WAIT_THREAD(f->thread_align);
          f->thread_align = 0;
        }
#endif

        /*
         * CREATE_THREAD(f->thread_align,thread_work_align,f);
         */

        /*
        *AlignAndFree(
        *        data_query, datasize_query,
        *        chunk_data_text, chunk_datasize_text,
        *        1,
        *        f);
        *
        *STATS_ADD_CLOCK(f, clock_align);
        */

        // try iterate through all tuples
        fprintf(OUTSTREAM, "%s\t%s\t%ld\t", gp_chunkname_query[f->j_chunk],
                gp_chunkname_text[f->i_chunk], f->reverse);

        tl_ = tl_prev_ ->next;


        while (tl_ != NULL) {
            fprintf(OUTSTREAM, "\t");
            t_ = tl_->first_tuple;
            int number_of_tuples = 0;
            while (t_ != NULL) {
                number_of_tuples++;
                fprintf(OUTSTREAM, "%ld %ld %ld,",
                        t_->occurrence,
                        t_->diagonal,
                        t_->leftsize);
                tl_prev_ = tl_;
                t_ = t_->next;
            }

            //nexttuplelist:
            tl_prev_ = tl_;
            tl_ = tl_->next;
        }
        fprintf(OUTSTREAM, "\n");
        //fprintf(OUTSTREAM, "%ld\n", number_of_tuples);

    }/*end of  [4] (chunk loop) */



    /* >> */
    FREE(last_chunk_diag_used, LAST_CHUNK_DIAG_NB * sizeof(long int));
    /* << */


#ifdef CACHE
    FREE(last_tuple_pos_with_diag_hash ,((1 + HASH_DIAG(MEMOPT_TABSIZE)) * sizeof(long int));
#endif




    FREE_TAB(last_tuple_pos_with_diag, sizeof(long int));
    FREE_TAB(last_tuple_ptr_with_diag, sizeof(tuple *));
    {
        long int seed;
        for (seed = 0; seed < gp_nb_seeds; seed++) {
            FREE(keylist_query[seed], keylist_query_size[seed]);
            FREE(first_query[seed], first_query_size[seed]);
        }
    }
    return 0;
}



