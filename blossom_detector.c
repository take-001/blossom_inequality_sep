#include <stdio.h>
#include <stdlib.h>
#include "machdefs.h"
#include "util.h"
#include "tsp.h"
#include "macrorus.h"

// Revised Blossom Loop Function
int revised_blossom_loop(int ncount, int ecount, int *elist, double *x, int silent, CCrandstate *rstate) {
    int rval = 0;
    int cutcount, cut_added;
    int inside = 0, outside = 0;
    int num_loop = 1; // Number of outer loops for separation routines
    CCtsp_lpcut_in *cuts = NULL;

    // Validate input
    if (ncount <= 0 || ecount <= 0 || !elist || !x) {
        fprintf(stderr, "Invalid input to revised_blossom_loop\n");
        return 1;
    }

    printf("Debug: Input Validation Passed\n");
    printf("Debug: ncount = %d, ecount = %d\n", ncount, ecount);



    do {
        do {
            cut_added = 0;

            // Fast Blossoms
            printf("Debug: Starting Fast Blossoms\n");
            rval = CCtsp_fastblossom(&cuts, &cutcount, ncount, ecount, elist, x);
            if (rval) {
                fprintf(stderr, "CCtsp_fastblossom failed\n");
                goto CLEANUP;
            }
            printf("Debug: Fast Blossoms completed with %d cuts\n", cutcount);

            if (cutcount) {
                cut_added += cutcount;
            }

            // Groetschel-Holland Fast Blossoms
            printf("Debug: Starting Groetschel-Holland Fast Blossoms\n");
            rval = CCtsp_ghfastblossom(&cuts, &cutcount, ncount, ecount, elist, x);
            if (rval) {
                fprintf(stderr, "CCtsp_ghfastblossom failed\n");
                goto CLEANUP;
            }
            printf("Debug: Groetschel-Holland Fast Blossoms completed with %d cuts\n", cutcount);

            if (cutcount) {
                cut_added += cutcount;
            }

            // Exact Blossoms
            printf("Debug: Starting Exact Blossoms\n");
            rval = CCtsp_exactblossom(&cuts, &cutcount, ncount, ecount, elist, x, rstate);
            if (rval) {
                fprintf(stderr, "CCtsp_exactblossom failed\n");
                goto CLEANUP;
            }
            printf("Debug: Exact Blossoms completed with %d cuts\n", cutcount);

            if (cutcount) {
                cut_added += cutcount;
            }

        } while (cut_added > 0);

        outside++;
    } while (outside < num_loop);

CLEANUP:
    while (cuts) {
        CCtsp_lpcut_in *next = cuts->next;
        CCtsp_free_lpcut_in(cuts);
        cuts = next;
    }
    return rval;
}

// Main Function
int main(int argc, char **argv) {
    int ncount = 500; // Example: 50 nodes
    int ecount = 1000; // Example: 100 edges

    // Dynamically allocate edge list and fractional values
    int *elist = malloc(2 * ecount * sizeof(int));
    double *x = malloc(ecount * sizeof(double));
    if (!elist || !x) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    // Generate example edge list and fractional values
    for (int i = 0; i < ecount; i++) {
        elist[2 * i] = rand() % ncount;           // Random node 1
        elist[2 * i + 1] = rand() % ncount;       // Random node 2
        if (elist[2 * i] == elist[2 * i + 1]) {   // Avoid self-loops
            elist[2 * i + 1] = (elist[2 * i + 1] + 1) % ncount;
        }
        x[i] = (rand() % 100) / 100.0; // Fractional value between 0 and 1
    }

    int silent = 0;
    CCrandstate rstate;
    CCutil_sprand(12345, &rstate); // Initialize random state

    printf("Running revised blossom loop...\n");
    int rval = revised_blossom_loop(ncount, ecount, elist, x, silent, &rstate);
    if (rval) {
        fprintf(stderr, "Blossom loop failed\n");
        free(elist);
        free(x);
        return 1;
    }
    printf("Blossom loop completed successfully.\n");

    free(elist);
    free(x);
    return 0;
}
