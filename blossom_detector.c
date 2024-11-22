#include <stdio.h>
#include <stdlib.h>
#include "machdefs.h"
#include "util.h"
#include "tsp.h"
#include "macrorus.h"

// Function Prototypes
int revised_blossom_loop(int ncount, int ecount, int *elist, double *x, int silent, CCrandstate *rstate);
void print_cuts(CCtsp_lpcut_in *cuts);
void free_cuts(CCtsp_lpcut_in *cuts);
void generate_fractional_solution(int ncount, int ecount, int *elist, double *x);

// Main Function
int main(int argc, char **argv) {
    int ncount = 10; // Number of nodes
    int ecount = ncount*(ncount)-1; // Number of edges for 2-regular graph

    // Dynamically allocate edge list and fractional values
    int *elist = malloc(2 * ecount * sizeof(int));
    double *x = malloc(ecount * sizeof(double));
    if (!elist || !x) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    // Generate edge list and fractional values ensuring constraints
    generate_fractional_solution(ncount, ecount, elist, x);

    // Print generated edges and fractional values
    printf("Generated Edge List and Fractional Values:\n");
    for (int i = 0; i < ecount; i++) {
        printf("Edge %d-%d: Fractional Value = %.4f\n", elist[2 * i], elist[2 * i + 1], x[i]);
    }

    int silent = 0;
    CCrandstate rstate;
    CCutil_sprand(12345, &rstate); // Initialize random state

    printf("\nRunning revised blossom loop...\n");
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

// Generate Edge List and Fractional Values
void generate_fractional_solution(int ncount, int ecount, int *elist, double *x) {
    double *vertex_sum = calloc(ncount, sizeof(double));
    if (!vertex_sum) {
        fprintf(stderr, "Memory allocation failed for vertex sum array.\n");
        exit(1);
    }

    // Generate edges ensuring all vertices are connected in a cycle
    for (int i = 0; i < ncount; i++) {
        elist[2 * i] = i;
        elist[2 * i + 1] = (i + 1) % ncount; // Connect to next vertex, wrapping around
    }

    // Assign fractional values ensuring sum of incident edges equals 2
    for (int i = 0; i < ecount; i++) {
        int u = elist[2 * i];
        int v = elist[2 * i + 1];

        // Calculate the remaining fractional capacity for each vertex
        double remaining_u = 2.0 - vertex_sum[u];
        double remaining_v = 2.0 - vertex_sum[v];

        // Determine the maximum assignable value to this edge
        double max_assignable = (remaining_u < remaining_v) ? remaining_u : remaining_v;
        max_assignable = (max_assignable > 1.0) ? 1.0 : max_assignable; // Cap at 1.0

        if (i == ecount - 1 || max_assignable <= 0.0) {
            // Final adjustment for the last edge or when max_assignable is 0
            x[i] = (remaining_u < remaining_v) ? remaining_u : remaining_v;
        } else {
            // Random fractional value between 0 and max_assignable
            x[i] = ((double)rand() / RAND_MAX) * max_assignable;
        }

        // Update the vertex sum for each incident vertex
        vertex_sum[u] += x[i];
        vertex_sum[v] += x[i];
    }

    // Debugging: Print the sum of fractional values for each vertex
    printf("Vertex fractional sums:\n");
    for (int i = 0; i < ncount; i++) {
        printf("Vertex %d: Sum = %.4f\n", i, vertex_sum[i]);
    }

    free(vertex_sum);
}



// Revised Blossom Loop
int revised_blossom_loop(int ncount, int ecount, int *elist, double *x, int silent, CCrandstate *rstate) {
    int rval = 0;
    int cutcount, cut_added;
    int outside = 0;
    int num_loop = 10;
    CCtsp_lpcut_in *cuts = NULL;

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
            print_cuts(cuts);
            free_cuts(cuts);
            cuts = NULL;

            // Groetschel-Holland Fast Blossoms
            printf("Debug: Starting Groetschel-Holland Fast Blossoms\n");
            rval = CCtsp_ghfastblossom(&cuts, &cutcount, ncount, ecount, elist, x);
            if (rval) {
                fprintf(stderr, "CCtsp_ghfastblossom failed\n");
                goto CLEANUP;
            }
            printf("Debug: Groetschel-Holland Fast Blossoms completed with %d cuts\n", cutcount);
            print_cuts(cuts);
            free_cuts(cuts);
            cuts = NULL;

            // Exact Blossoms
            printf("Debug: Starting Exact Blossoms\n");
            rval = CCtsp_exactblossom(&cuts, &cutcount, ncount, ecount, elist, x, rstate);
            if (rval) {
                fprintf(stderr, "CCtsp_exactblossom failed\n");
                goto CLEANUP;
            }
            printf("Debug: Exact Blossoms completed with %d cuts\n", cutcount);
            print_cuts(cuts);
            free_cuts(cuts);
            cuts = NULL;

        } while (cut_added > 0);

        outside++;
    } while (outside < num_loop);

CLEANUP:
    if (cuts) {
        free_cuts(cuts);
        cuts = NULL;
    }
    return rval;
}

// Print Cuts
void print_cuts(CCtsp_lpcut_in *cuts) {
    CCtsp_lpcut_in *current_cut = cuts;
    int comb_index = 0;

    while (current_cut) {
        printf("Comb %d:\n", ++comb_index);

        // Iterate through all cliques in the comb
        for (int t = 0; t < current_cut->cliquecount; t++) {
            for (int i = 0; i < current_cut->cliques[t].segcount; i++) {
                printf("[%d, %d] ", 
                       current_cut->cliques[t].nodes[i].lo, 
                       current_cut->cliques[t].nodes[i].hi);
            }
            printf("\n");
        }

        current_cut = current_cut->next;
    }
}


// Free Cuts
void free_cuts(CCtsp_lpcut_in *cuts) {
    while (cuts) {
        CCtsp_lpcut_in *next = cuts->next;
        CCtsp_free_lpcut_in(cuts);
        cuts = next;
    }
    printf("freeing cuts completed\n");
}
