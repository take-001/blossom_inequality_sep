#include "machdefs.h"
#include "util.h"
#include "tsp.h"
#include "macrorus.h"

// Revised Blossom Loop Function
int revised_blossom_loop(int ncount, int ecount, int *elist, double *x, int silent, CCrandstate *rstate) {
    int rval = 0;
    int cutcount, cut_added;
    int inside = 0, outside = 0;
    double z;
    CCtsp_lpcut_in *cuts = NULL;
    int num_loop = 1;

    // Validate input
    if (ncount <= 0 || ecount <= 0 || !elist || !x) {
        fprintf(stderr, "Invalid input to revised_blossom_loop\n");
        return 1;
    }

    do {
        do {
            cut_added = 0;

            // Fast Blossoms
            printf("Debug: Starting Fast Blossoms\n");
            // CCutil_start_timer(NULL); // Commented out
            rval = CCtsp_fastblossom(&cuts, &cutcount, ncount, ecount, elist, x);
            if (rval) {
                fprintf(stderr, "CCtsp_fastblossom failed\n");
                goto CLEANUP;
            }
            // z = CCutil_stop_timer(NULL, 0); // Commented out
            printf("Debug: Fast Blossoms completed\n");

            if (!silent) {
                printf("Found %d Fast Blossoms\n", cutcount);
            }
            if (cutcount) {
                cut_added += cutcount;
            }

            // Groetschel-Holland Fast Blossoms
            printf("Debug: Starting Groetschel-Holland Fast Blossoms\n");
            // CCutil_start_timer(NULL); // Commented out
            rval = CCtsp_ghfastblossom(&cuts, &cutcount, ncount, ecount, elist, x);
            if (rval) {
                fprintf(stderr, "CCtsp_ghfastblossom failed\n");
                goto CLEANUP;
            }
            // z = CCutil_stop_timer(NULL, 0); // Commented out
            printf("Debug: Groetschel-Holland Fast Blossoms completed\n");

            if (!silent) {
                printf("Found %d Groetschel-Holland Blossoms\n", cutcount);
            }
            if (cutcount) {
                cut_added += cutcount;
            }

            // Exact Blossoms
            printf("Debug: Starting Exact Blossoms\n");
            // CCutil_start_timer(NULL); // Commented out
            rval = CCtsp_exactblossom(&cuts, &cutcount, ncount, ecount, elist, x, rstate);
            if (rval) {
                fprintf(stderr, "CCtsp_exactblossom failed\n");
                goto CLEANUP;
            }
            // z = CCutil_stop_timer(NULL, 0); // Commented out
            printf("Debug: Exact Blossoms completed\n");

            if (!silent) {
                printf("Found %d Exact Blossoms\n", cutcount);
            }
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
    int ncount = 5; // Example: 5 nodes
    int ecount = 8; // Example: 8 edges
    int elist[] = {0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3, 2, 4, 3, 4}; // Example edges
    double x[] = {0.5, 0.3, 0.7, 0.2, 0.1, 0.6, 0.4, 0.8}; // Example fractional values
    int silent = 0;
    CCrandstate rstate;

    CCutil_sprand(12345, &rstate); // Initialize random state

    printf("Running revised blossom loop...\n");
    int rval = revised_blossom_loop(ncount, ecount, elist, x, silent, &rstate);
    if (rval) {
        fprintf(stderr, "Blossom loop failed\n");
        return 1;
    }
    printf("Blossom loop completed successfully.\n");

    return 0;
}
