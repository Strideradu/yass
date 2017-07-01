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
#include <errno.h>
#include <math.h>

/* 1) include utils macro */
#include "util.h"
/* 2) include global variables */
#include "global_var.h"
/* 3) the current file defs */
#include "kword.h"




/*
 *   Build a  DNA bases table  (code 0..3) from a DNA file
 *
 * - Input  : "filename"                                    : DNA file name
 * - Output : "p_data"                                      : return the allocated data table pointer
 *            "p_datasize"                                  : return the allocated data table size
 *            "p_nbchunks","p_chunkname", and "p_chunksize" : return respectively the number of chunks fond in the multifasta file, their names and their size
 *                                                            (set to NULL if you do not need these informations
 */



#define WORDCOUNTANDSTAT(count, letter, word) {(word) <<= 2;  (word) |= (letter); (word) &= 0xfc00003f;  if (((word) & 0x3f) == (word)) count[(word)]++;}


long int CreateData(/* in */  char *filename,
        /* out */ char **p_data,         /* out */ long int *p_datasize,
        /* out */ long int *p_nbchunks,  /* out */ char ***p_chunkname,
        /* out */ long int **p_chunksize,/* out */ long int **p_chunkstart,
        /* out */ long int *nb_letters  /* [4]  */ ,
        /* out */ long int *nb_triplets /* [64] */ ) {

    long int i = 0;
    long int j = 0;
    int car = 0;
    long int datasize = 0;
    long int cumulength = 0;
    char *data = NULL;
    FILE *file = NULL;

    long int nbchunks = 0;
    char **chunkname = NULL;
    long int *chunksize = NULL;
    long int *chunkstart = NULL;
    long unsigned int word = 0xfc000000;

    long int line = 1, column = 1;
    int invalid_report[256] = {0};
    char chunkname_buffer[CHUNKNAMESIZE + 1];

    /* (1) count number of DNA bases availaible */


/*
 *      A --> adenosine           K --> G T (keto)
 *      C --> cytidine            M --> A C (amino)
 *      G --> guanine             S --> G C (strong)
 *      T --> thymidine           W --> A T (weak)
 *      R --> G A (purine)        V --> G C A
 *      Y --> T C (pyrimidine)    B --> G T C
 *      N --> A G C T (any)       D --> G A T
 *      U --> uridine             H --> A C T
 */

    if (!(file = fopen(filename, "r"))) {
        fprintf(stderr, "* Error : unable to open \"%s\": errno = %d \n", filename, errno);
        exit(1);
    }

    /* comment */
    while ((car = fgetc(file)) != EOF) {
        switch (car) {
            case 'A':
            case 'a':
            case 'C':
            case 'c':
            case 'G':
            case 'g':
            case 'T':
            case 't':
            case 'R':
            case 'r':
            case 'Y':
            case 'y':
            case 'N':
            case 'n':
            case 'U':
            case 'u':
            case 'K':
            case 'k':
            case 'M':
            case 'm':
            case 'S':
            case 's':
            case 'W':
            case 'w':
            case 'V':
            case 'v':
            case 'B':
            case 'b':
            case 'D':
            case 'd':
            case 'H':
            case 'h':
                column++;
                datasize++;
                break;
            case '>':
            case '#':
            case ':':
                if (column == 1) {
                    column++;
                    nbchunks++;
                    do {
                        /* seek eol */
                        car = fgetc(file);
                        column++;
                    } while (car != EOF && car != '\n' && car != '\r');
                    column = 1;
                    line++;
                } else {
                    fprintf(stderr, "* Error : file \"%s\" : invalid header char \'%c\' on line %ld, column %ld\n",
                            filename, car, line, column);
                    exit(1);
                }
                break;
            case '\r':
            case '\n':
                column = 1;
                line++;
                break;
            case ' ':
            case '\t':
                column++;
                break;
            default:
                if (!invalid_report[(unsigned) car]) {
                    invalid_report[(unsigned) car] = 1;
                    fprintf(stderr,
                            "* Warning : file \"%s\" : invalid character \'%c\' (\"0x%2X\") on line %ld, column %ld (ignored and not reported anymore on this file)\n",
                            filename, car, (unsigned) car, line, column);
                }
                column++;
                break;
        } /*switch */
    } /* while */
    if (datasize < 0) {
        fprintf(stderr, "* Error : file \"%s\" is too large\n", filename);
        exit(1);
    }

    if (!datasize) {
        fprintf(stderr, "* Error : empty file \"%s\"\n", filename);
        exit(1);
    }

    /* datasize known */
    fseek(file, 0L, SEEK_SET);

    /* (2) allocate one memory blok for the full DNA data */
    data = (char *) MALLOC(datasize * sizeof(char));
    ASSERT(data, CreateData);

    if (nbchunks > 0) {
        chunkname = (char **) MALLOC(nbchunks * sizeof(char *));
        ASSERT(chunkname, CreateData);
        memset(chunkname, 0, nbchunks * sizeof(char *));
        chunksize = (long int *) MALLOC(nbchunks * sizeof(long int));
        ASSERT(chunksize, CreateData);
        memset(chunksize, 0, nbchunks * sizeof(long int));
        chunkstart = (long int *) MALLOC(nbchunks * sizeof(long int));
        ASSERT(chunkstart, CreateData);
        memset(chunkstart, 0, nbchunks * sizeof(long int));
    }


    /* (3) fill the DNA data */
    /* description of sentance */
    while ((car = fgetc(file)) != EOF) {

        switch (car) {

            case 'A':
                data[i++] = (char) 0;
                nb_letters[0]++;
                WORDCOUNTANDSTAT(nb_triplets, 0, word);
                break;

            case 'C':
                data[i++] = (char) 1;
                nb_letters[1]++;
                WORDCOUNTANDSTAT(nb_triplets, 1, word);
                break;

            case 'G':
                data[i++] = (char) 2;
                nb_letters[2]++;
                WORDCOUNTANDSTAT(nb_triplets, 2, word);
                break;

            case 'T':
                data[i++] = (char) 3;
                nb_letters[3]++;
                WORDCOUNTANDSTAT(nb_triplets, 3, word);
                break;

            case 'R':
                data[i++] = (char) 4;
                break;

            case 'Y':
                data[i++] = (char) 5;
                break;

            case 'N':
                data[i++] = (char) 6;
                break;

            case 'U':
                data[i++] = (char) 7;
                nb_letters[3]++;
                WORDCOUNTANDSTAT(nb_triplets, 3, word);
                break;

            case 'K':
                data[i++] = (char) 8;
                break;

            case 'W':
                data[i++] = (char) 9;
                break;

            case 'S':
                data[i++] = (char) 10;
                break;

            case 'M':
                data[i++] = (char) 11;
                break;

            case 'V':
                data[i++] = (char) 12;
                break;

            case 'H':
                data[i++] = (char) 13;
                break;

            case 'D':
                data[i++] = (char) 14;
                break;

            case 'B':
                data[i++] = (char) 15;
                break;

            case 'a':
                data[i++] = (char) 16;
                nb_letters[0]++;
                WORDCOUNTANDSTAT(nb_triplets, 0, word);
                break;

            case 'c':
                data[i++] = (char) 17;
                nb_letters[1]++;
                WORDCOUNTANDSTAT(nb_triplets, 1, word);
                break;

            case 'g':
                data[i++] = (char) 18;
                nb_letters[2]++;
                WORDCOUNTANDSTAT(nb_triplets, 2, word);
                break;

            case 't':
                data[i++] = (char) 19;
                nb_letters[3]++;
                WORDCOUNTANDSTAT(nb_triplets, 3, word);
                break;

            case 'r':
                data[i++] = (char) 20;
                break;

            case 'y':
                data[i++] = (char) 21;
                break;

            case 'n':
                data[i++] = (char) 22;
                break;

            case 'u':
                data[i++] = (char) 23;
                nb_letters[3]++;
                WORDCOUNTANDSTAT(nb_triplets, 3, word);
                break;

            case 'k':
                data[i++] = (char) 24;
                break;

            case 'w':
                data[i++] = (char) 25;
                break;

            case 's':
                data[i++] = (char) 26;
                break;

            case 'm':
                data[i++] = (char) 27;
                break;

            case 'v':
                data[i++] = (char) 28;
                break;

            case 'h':
                data[i++] = (char) 29;
                break;

            case 'd':
                data[i++] = (char) 30;
                break;

            case 'b':
                data[i++] = (char) 31;
                break;

            case '\n':
            case '\r':
            case ' ':
            case '\t':
                break;
            case '>':
            case '#':
            case ':':
                /* --[A] chunk name / size update */
                if (j > 0) {
                    /* update previous length if needed */
                    if (i > cumulength)
                        chunksize[j - 1] = i - cumulength;
                    else
                        /* otherwise, no data between : rewrite the new chunk */
                        j--;
                } else {
                    /* j == 0  this is the first chunk */
                    /* "incorrect file format" : remove the "dust" before the first chunk */
                    i = 0;
                }
                word = 0xfc000000;
                cumulength = i;
                chunkstart[j] = i;
                {
                    int c = 0;
                    do {
                        /* seek eol/eof */
                        car = fgetc(file);
                        if (car == ' ' || car == '\t')
                            car = '_';
                        if (car != EOF && car != '\n' && car != '\r' && c < CHUNKNAMESIZE)
                            chunkname_buffer[c++] = car;
                    } while (car != EOF && car != '\n' && car != '\r');
                    chunkname_buffer[c] = '\0';
                    chunkname[j] = (char *) MALLOC((c + 1) * sizeof(char));
                    ASSERT(chunkname[j], CreateData);
                    strncpy(chunkname[j], chunkname_buffer, c);/*FIXME*/
                    chunkname[j][c] = '\0';
                }
                j++;
                break;
            default:
                /* ignore */;
        } /* switch */
    } /* while */

    fclose(file);

    /* --[B] last chunk size updated */
    if (j > 0) {
        if (i > cumulength)
            chunksize[j - 1] = i - cumulength;
        else
            /* otherwise, no data between : rewrite the new chunk */
            j--;
    }

    nbchunks = j;

    /* --[C] Raw format : simulate one big fasta chunk */
    if (nbchunks <= 0) {
        nbchunks = 1;
        chunkname = (char **) MALLOC(nbchunks * sizeof(char *));
        ASSERT(chunkname, CreateData);
        chunksize = (long int *) MALLOC(nbchunks * sizeof(long int));
        ASSERT(chunksize, CreateData);
        chunkstart = (long int *) MALLOC(nbchunks * sizeof(long int));
        ASSERT(chunkstart, CreateData);
        chunkname[0] = filename;
        chunksize[0] = datasize;
        chunkstart[0] = 0;
    }
#ifdef DEBUG_CREATEDATA
    for (i = 0; i < nbchunks; i++)
      fprintf(stderr,
              " chunksize = %ld chunkstart = %ld chunkname = \"%s\" \n",
              chunksize[i], chunkstart[i], chunkname[i]);
#endif

    /* (4) */
    /* FIXME : set some low complexity filter here ... */



    /* (5) return allocated and computed data */
    *p_datasize = datasize;
    *p_data = data;
    if (p_nbchunks != NULL)
        *p_nbchunks = nbchunks;
    if (p_chunkname != NULL)
        *p_chunkname = chunkname;
    if (p_chunksize != NULL)
        *p_chunksize = chunksize;
    if (p_chunkstart != NULL)
        *p_chunkstart = chunkstart;
    return 0;
}


/*
 *   Build the complementary table of a DNA mfasta file
 *
 * - Input  : "data"                                        : DNA data
 *            "datasize"                                    : DNA data size (full length)
 *            "nbchunks"                                    : number of blocks inside the mfasta file
 *            "chunksize[nbchunks]"                         : size of each block (sum must be equal to "datasize")
 *            "chunkstart[nbchunks]"                        : start of each block (shift according to the begining)
 * - Output : "p_data_rev_out"                              : return the allocated data table pointer
 *                                                            (each block has been reversed, but keep the same chunkstart "position")
 */


long int CreateReverseComplement(
        /* in  */ char *data,
        /* in  */ long int datasize,
        /* in  */ long int nbchunks,
        /* in  */ long int *chunksize,
        /* in  */ long int *chunkstart,
        /* out */ char **p_data_rev_out) {
    long int i;
    char *data_rev = (char *) MALLOC(datasize * sizeof(char));
    ASSERT(data_rev, CreateDataReverseComplement);
    for (i = 0; i < nbchunks; i++) {
        long int currentchunksize = chunksize[i];
        char *a = data + chunkstart[i];
        char *b = data_rev + chunkstart[i] + currentchunksize;
        long int i = 0;
        for (i = 0; i < currentchunksize; i++)
            *(--b) = (char) COMPLEMENT(*a++);
    }
    *p_data_rev_out = data_rev;
    return 0;
}


/*
 *   Build the table that give a number of dots to each Qchunk
 *
 * - Input  : "nbchunks"                                    : number of blocks inside the mfasta file
 *            "chunksize[nbchunks]"                         : size of each block
 * - Return : a table giving how much dots have to be printed ...
 */

long int *ComputeDotsTable(long int nbchunks, long int *chunksize) {

    long int cumulativesize = 0;
    long int overallsize = 0;
    long int nbdots = 0;
    long int *dots = (long int *) MALLOC(nbchunks * sizeof(long int));
    memset(dots, '\0', nbchunks * sizeof(long int));
    ASSERT(dots, ComputeDotsTable);

    {
        long int i;
        for (i = 0; i < nbchunks; i++)
            overallsize += chunksize[i];
    }

    {
        long int i;
        for (i = 0; i < nbchunks; i++) {
            cumulativesize += chunksize[i];
            while (nbdots * overallsize < ((gp_reverse == 2) ? (NBL / 2) : (NBL)) * cumulativesize) {
                nbdots++;
                dots[i]++;
            }
        }
    }
    return dots;
}


/*
 * seed pattern size computation
 * (gp_seed_bitweight, gp_seed_span)
 */


long int ComputeLengthAndSortSeeds() {
    long int seed;
    long int i;
    gp_seeds_span_max = 0;
    gp_seeds_span_min = INFINITY_INT;

    /* 1) compute length */
    for (seed = 0; seed < gp_nb_seeds; seed++) {
        gp_seeds_bitweight[seed] = 0;
        gp_seeds_span[seed] = 0;

        for (i = 0; i < (long int) strlen(gp_motifs[seed]); i++) {
            if (gp_motifs[seed][i] == '#') {
                gp_seeds_bitweight[seed] += 2;
            } else {
                if (gp_motifs[seed][i] == '@') {
                    gp_seeds_bitweight[seed] += 1;
                } else {
                    /* if another char is given : replace it by - */
                    gp_motifs[seed][i] = '-';
                }
            }
            gp_seeds_span[seed]++;
        }

        gp_motifs[seed][gp_seeds_span[seed]] = '\0';
        gp_seeds_span_max = MAX(gp_seeds_span_max, gp_seeds_span[seed]);
        gp_seeds_span_min = MIN(gp_seeds_span_min, gp_seeds_span[seed]);
    }

    /* 2) sort seeds by increasing span (thus the maximal span will always match after the others) */
    {
        long int i;
        for (i = 1; i < gp_nb_seeds; i++) {
            long int j;
            for (j = 0; j < i; j++) {
                if (gp_seeds_span[j] > gp_seeds_span[i]) {
                    /* swap */
                    long int s = gp_seeds_span[i];
                    long int w = gp_seeds_bitweight[i];
                    char *m = gp_motifs[i];
                    gp_seeds_span[i] = gp_seeds_span[j];
                    gp_seeds_bitweight[i] = gp_seeds_bitweight[j];
                    gp_motifs[i] = gp_motifs[j];
                    gp_seeds_span[j] = s;
                    gp_seeds_bitweight[j] = w;
                    gp_motifs[j] = m;

                }
            }
        }
    }
    return 0;
}


/*
 * remove from the "first" index
 * seed codes that are low complexity ones
 * (only one letter, or two letters)
 */


static void BuildRedMask(long int seed, long int *first) {
    long int mask_copy = 0;
    long int c_log = 0;
    /* precompute some seed positions */
    {
        long int i;
        for (i = 0; i < gp_seeds_span[seed]; i++) {
            char seed_element = gp_motifs[seed][i];
            switch (seed_element) {
                case '#' :
                    mask_copy <<= 1;
                    mask_copy |= 1;
                    c_log++;
                    break;
                case '@' :
                    mask_copy <<= 1;
                    c_log++;
                    break;
            }
        }
    }

    /* (weight >= 6) => for all single letters */
    if (gp_seeds_bitweight[seed] >= 6 * 2) {
        long int a;
        for (a = 0; a < 4; a++) {

            /* build all words containing only 'a' letters */
            long int keycode = 0;
            long int i = 0;

            /* build (a)-key code */
            for (i = 0; i < c_log; i++) {
                long int i_pow = 1 << i;
                if (i_pow & mask_copy) {
                    keycode <<= 2;
                    keycode |= a;
                } else {
                    keycode <<= 1;
                    keycode |= a & 1;
                }
            }
            /* reset this keycode */
            first[keycode] = -2;
        }
    }

    /* (weight >= 8) => for all pair of letters */
    if (gp_seeds_bitweight[seed] >= 8 * 2) {
        long int a;
        for (a = 1; a < 4; a++) {
            long int b;
            for (b = 0; b < a; b++) {
                long int c;

                /* build all words containing only 'a' and 'b' letters */
                for (c = 0; c < (1 << c_log); c++) {
                    long int keycode = 0;
                    long int i = 0;

                    /* build (a,b)-key code */
                    for (i = 0; i < c_log; i++) {
                        long int i_pow = 1 << i;
                        if (i_pow & mask_copy) {
                            keycode <<= 2;
                            keycode |= ((c & i_pow) ? a : b);
                        } else {
                            keycode <<= 1;
                            keycode |= ((c & i_pow) ? a : b) & 1;
                        }
                    }

                    /* reset this keycode */
                    first[keycode] = -2;
                }
            }
        }
    }
}


/*
 * Create a linkedlist of kwords with param "list" and index "first"
 *
 * - Input : "data"    : dna input table
 *           "datasize": dna input table size
 *           "seed"    : seed number
 * - Output: "p_first" : return a table of size (2^w). For each code w, it contains the first
 *                       position occurence of the seed with such code on the sequence
 *           "p_list"  : return a table of size (keysize). For each position i, it contains next
 *                       position j such that seed i and j code are equal.
 */


#define ISTART 0

long int CreateKeyList( /* in */  char *data,         /* in  */ long int datasize,      /* in */  long int seed,
        /* out */ long int **p_first, /* out */ long int *p_first_size,
        /* out */ long int **p_list,  /* out */ long int *p_list_size, /* out */ Feature *feature) {


    long int *bucklet = NULL;
    long int *first = NULL;
    long int *list = NULL;
    long int *count = NULL;
#ifdef KEYLISTCOMPRESS
    long int *listcompress = NULL;
#endif
    long int sizeofbucklet = 1;

    /* prune algorithm */
    long int maxchainsize = 0;
    long int nbkeys = 0;
    long int fastcountall = 0;
    long int fastcounter[sizeof(unsigned long int) * 8] = {0};
    long int imask[sizeof(unsigned long int) * 8] = {0};

    /* (1) allocation of 3 table */
    sizeofbucklet <<= (gp_seeds_bitweight[seed]);

    bucklet = (long int *) MALLOC(sizeofbucklet * sizeof(long int));
    ASSERT(bucklet, CreateKeyList);
    memset(bucklet, 0x00, sizeofbucklet * sizeof(long int));

    first = (long int *) MALLOC((sizeofbucklet + 1) * sizeof(long int));
    ASSERT(first, CreateKeyList);
    memset(first, 0xff, sizeofbucklet * sizeof(long int));

    count = (long int *) MALLOC(sizeofbucklet * sizeof(long int));
    ASSERT(count, CreateKeyList);
    memset(count, 0x0, sizeofbucklet * sizeof(long int));

    list = (long int *) MALLOC((datasize + 1) * sizeof(long int));
    ASSERT(list, CreateKeyList);
    memset(list, 0xff, (datasize + 1) * sizeof(long int));


    /* (2) list creation algorithm */
    BuildRedMask(seed, first); /*
                             * redmask (-2) is for unallowed key = "0, all G, all C, all A, all T"
                             * and dinucleotides only composed keys ....
                             */
    {
        long int pos;
        for (pos = 0; pos <= datasize - gp_seeds_span[seed]; pos++) {
            long int pos_end, w;
            KEY(w, pos_end, data, pos, seed);
            if (w >= 0) {
                if (first[w] != -2) {
                    if (first[w] == -1) {
                        first[w] = pos_end;
                    } else {
                        list[bucklet[w]] = pos_end;
                    }
                    bucklet[w] = pos_end;
                    count[w]++;
                    nbkeys++;
                    maxchainsize = MAX(count[w], maxchainsize);
                }
            }
        }
    }
    FREE(bucklet, sizeofbucklet * sizeof(long int));



    /* (3) median-mean selection */

    /* (3.0) init counters */
    imask[0] = 1;
    {
        long int i;
        for (i = 1; i < sizeof(long int) * 8; i++) {
            imask[i] = ((imask[i - 1]) << 1) | 1;
            fastcounter[i] = 0;
        }
    }

    /* (3.1) fast approximate counting */
    {
        long int w;
        for (w = 0; w < sizeofbucklet; w++) {
            long int i;
            for (i = ISTART; i < sizeof(long int) * 8; i++) {
                if ((imask[i] & count[w]) == count[w]) {
                    fastcounter[i]++;
                    break;
                }
            }
        }
    }

#ifdef DEBUGFASTCOUNT
    {
      long int i;
      for  (i = ISTART ; i < sizeof(long int)*8 ; i++ ) {
        fprintf(stderr, " fastcount[%ld] = %ld \n", i , fastcounter[i]);
      }
    }
#endif

    /* (3.2) select part to be removed part */
    {
        long int i;
        long int w;
        for (i = ISTART; i < sizeof(long int) * 8; i++) {
            fastcountall += fastcounter[i];
            if (fastcountall * 1.0 > sizeofbucklet * (1.0 - (gp_t + 1e-32) / sqrt((double) sizeofbucklet)))
                break;
        }

        /* (3.3) remove gp_t percent */


        for (w = 0; w < sizeofbucklet; w++) {
            /* normalize */
            if (first[w] < 0) {
                first[w] = -1;
            } else {
                if ((count[w] & imask[i]) != count[w]) {
                    first[w] = -1;
                    nbkeys -= count[w];
                    STATS_NB_KEYS_REMOVED_INC(feature);
#ifdef DEBUGRMKEY
                    fprintf(stderr,"- removing key : w=%lx, nbocc=%ld\n",w,count[w]);
#endif
                }
            }
        }
    }
    FREE(count, sizeofbucklet * sizeof(long int));

    /* (4) keylist compression */
#ifdef KEYLISTCOMPRESS
    listcompress = (long int *) MALLOC((nbkeys + 1) * sizeof(long int));
    ASSERT(listcompress, CreateKeyList);

    {
        long int i = 0;
        long int w;
        for (w = 0; w < sizeofbucklet; w++) {
            long int j = first[w];
            first[w] = i;
            while (j >= 0) {
                listcompress[i++] = j;
                j = list[j];
            }
        }
        first[sizeofbucklet] = i;
    }
#endif

    /* (5) list "copy" */
#ifdef KEYLISTCOMPRESS
    FREE(list, (datasize + 1) * sizeof(long int));
    *p_list = listcompress;
    *p_list_size = (nbkeys + 1) * sizeof(long int);
#else
    *p_list       = list;
    *p_list_size  = (datasize+1) * sizeof(long int);
#endif

    *p_first = first;
    *p_first_size = (sizeofbucklet + 1) * sizeof(long int);
    return 0;

}


/*
 * COUNTS
 */

long int lowercaseCount(char *data, long int pos, long int len) {
    long int i;
    long int c = 0;
    for (i = 0; i < len; i++)
        if (data[pos + i] >= 16)
            c++;
    return c;
}

long int nCount(char *data, long int pos, long int len) {
    long int i;
    long int c = 0;
    for (i = 0; i < len; i++)
        if ((data[pos + i] == 6) || (data[pos + i] == (16 + 6)))
            c++;
    return c;
}


/*
 * COMPLEXITY CHECK
 */

double Entropy(char *data, long int pos, long int len, long int k) {
    long int p4k = 1 << (k << 1);       /*4^k */
    long int mask = 0;
    long int key = 0;
    long int i = 0;
    double result = 0.0;
    long int *freq = (long int *) MALLOC(p4k * sizeof(long int));
    ASSERT(freq, Entropy);
    memset(freq, '\0', sizeof(long int) * p4k);

    for (i = 0; i < k; i++) {
        mask <<= 2;
        mask |= 0x3;
        key <<= 2;
        key |= backcode[(long int) data[pos + i]];
    }


    for (i = k; i < len; i++) {
        key <<= 2;
        key |= backcode[(long int) data[pos + i]];
        key &= mask;
        freq[key]++;
    }


    for (i = 0; i < p4k; i++) {
        if (freq[i] != 0) {
            result -=
                    ((double) (freq[i]) / (double) (len - k + 1)) *
                    log((double) (freq[i]) / (double) (len - k + 1));
        }
    }
    result /= log(2.0000);
    FREE(freq, p4k * sizeof(long int));
    return result;
}


/*
 * Display sequence fragments
 */

long int DisplaySequence(char *data, long int pos, long int len) {

    long int j;
    for (j = pos; j < pos + len; j++) {
        fprintf(stdout, "%c", lookup[(long int) data[j]]);
    }
    fprintf(stdout, "\n");
    return 0;
}


long int
DisplayData(char *filename, char *data, long int size, long int nbchunks, char **chunkname, long int *chunksize,
            long int *chunkstrt) {
    long int i;
    fprintf(stderr, "filename : \"%s\"\n", filename);
    fprintf(stderr, "- size : %ld\n", size);
    fprintf(stderr, "- nbchunks : %ld\n", nbchunks);
    for (i = 0; i < nbchunks; i++) {
        fprintf(stdout, ">\"%s\"|[start:%ld,size:%ld]\n", chunkname[i], chunkstrt[i], chunksize[i]);
        DisplaySequence(data, chunkstrt[i], chunksize[i]);
    }
    return 0;
}

long int DisplaySeeds() {
    long int i = 0;
    for (i = 0; i < gp_nb_seeds; i++) {
        fprintf(stderr, "- gp_motif[%ld] = \"%s\", span = %3ld , bitweight = %3ld\n", i, gp_motifs[i], gp_seeds_span[i],
                gp_seeds_bitweight[i]);
    }
    return 0;
}
