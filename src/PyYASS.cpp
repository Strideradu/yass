//
// Created by Admin on 2017/7/24.
//

//#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "util.h"
#include "global_var.h"
#include "kword.h"
#include "tuple.h"
#include "assemble.h"
#include "prdyn.h"
#include "proba.h"
#include "display.h"
#include "regroup.h"
#include "threads.h"

/*
 * fillmatrixmttrtv
 */

void fillmatrixmttrtv() {
    long int x, y;
    for (x = 0; x < 32; x++) {
        for (y = 0; y < 32; y++) {
            if (((x >= 0 && x < 6) || (x >= 16 && x < 22)) && ((y >= 0 && y < 6) || (y >= 16 && y < 22)))
                if ((x % 16) == (y % 16))
                    gp_substitution_matrix[x][y] = (int) gp_costs[0]; /* match */
                else if ((x % 2) == (y % 2))
                    gp_substitution_matrix[x][y] = (int) gp_costs[2]; /* transition */
                else
                    gp_substitution_matrix[x][y] = (int) gp_costs[1]; /* transversion */
            else
                gp_substitution_matrix[x][y] = (int) gp_costs[3]; /* other */

            gp_cost_max_substitution_matrix = MAX(
                    ABS(gp_substitution_matrix[x][y]),
                    gp_cost_max_substitution_matrix
            );
            /* symbols U (7) ~ T (3) */
            if (y == 7)
                gp_substitution_matrix[x][y] = gp_substitution_matrix[x][3];
            if (x == 7)
                gp_substitution_matrix[x][y] = gp_substitution_matrix[3][y];
            if (y == 7 && x == 7)
                gp_substitution_matrix[x][y] = gp_substitution_matrix[3][3];
            if (y == 23)
                gp_substitution_matrix[x][y] = gp_substitution_matrix[x][3];
            if (x == 23)
                gp_substitution_matrix[x][y] = gp_substitution_matrix[3][y];
            if (y == 23 && x == 23)
                gp_substitution_matrix[x][y] = gp_substitution_matrix[3][3];
        }
    }
}
/*
int find_hits(int error_percent, int gap_percent ){
    // init size 9 seed
    gp_motifs[0] = (char *) strdup("###-#@-##@##");
    gp_seeds_bitweight[0] = 18;
    gp_seeds_span_min = 9;
    gp_seeds_span_max = 9;
    gp_nb_seeds = 1;

    // Init scoring matrix
    gp_substitution_matrix = lint_directtable(32, 32);
    fillmatrixmttrtv();

    // Init semaphores
    INIT_MUTEX();INIT_QUERY_MUTEX();

    // [2] initialize parameters and statistics
    gp_mutations_percent = error_percent;
    gp_indels_percent = gap_percent;
    gp_mutations = (double) gp_mutations_percent / 100.00;
    gp_indels = (double) gp_indels_percent / 100.00;
    gp_alpha = (double) gp_alpha_percent / 100.00;
    gp_rho_stat = statistical_bound_of_waiting_time2(
            (1 - gp_mutations),
            gp_seeds_bitweight[0] / 2,
            gp_alpha
    );
    gp_delta_stat = statistical_bound_of_randomwalk2(gp_indels,
                                                     gp_rho_stat,
                                                     gp_alpha);

    gp_border = gp_rho_stat + gp_delta_stat + 1;

    // [3] chaining algorithm memory required
    initialise_deltashift();

}
 */

