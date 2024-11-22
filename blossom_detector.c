#include <stdio.h>
#include <stdlib.h>
#include "machdefs.h"
#include "util.h"
#include "tsp.h"
#include "macrorus.h"

// Function Prototypes
int revised_blossom_loop(int ncount, int ecount, int *elist, double *x, int silent, CCrandstate *rstate);
void verify_and_print_comb(CCtsp_lpcut_in *cuts, int ncount, int ecount, int *elist, double *x);
void free_cuts(CCtsp_lpcut_in *cuts);
void generate_fractional_solution(int ncount, int ecount, int *elist, double *x);

// Main Function
int main(int argc, char **argv) {
    int ncount = 400; // Number of nodes
    int ecount = ncount * (ncount - 1) / 2; // Number of edges for a complete graph

    // Dynamically allocate edge list and fractional values
    int *elist = malloc(2 * ecount * sizeof(int));
    double *x = malloc(ecount * sizeof(double));
    if (!elist || !x) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    // Generate edge list and fractional values ensuring constraints
    generate_fractional_solution(ncount, ecount, elist, x);

    printf("\nRunning revised blossom loop...\n");
    int silent = 0;
    CCrandstate rstate;
    CCutil_sprand(12345, &rstate); // Initialize random state

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

// Generate Fractional Solution
void generate_fractional_solution(int ncount, int ecount, int *elist, double *x) {
    double *vertex_sum = calloc(ncount, sizeof(double));
    if (!vertex_sum) {
        fprintf(stderr, "Memory allocation failed for vertex sum array.\n");
        exit(1);
    }

    // Generate a complete graph
    int edge_index = 0;
    for (int i = 0; i < ncount; i++) {
        for (int j = i + 1; j < ncount; j++) {
            elist[2 * edge_index] = i;
            elist[2 * edge_index + 1] = j;
            edge_index++;
        }
    }

    // Assign fractional values
    for (int i = 0; i < ecount; i++) {
        int u = elist[2 * i];
        int v = elist[2 * i + 1];
        double remaining_u = 2.0 - vertex_sum[u];
        double remaining_v = 2.0 - vertex_sum[v];
        double max_assignable = (remaining_u < remaining_v) ? remaining_u : remaining_v;
        max_assignable = (max_assignable > 1.0) ? 1.0 : max_assignable;

        x[i] = ((double)rand() / RAND_MAX) * max_assignable;
        vertex_sum[u] += x[i];
        vertex_sum[v] += x[i];
    }


    free(vertex_sum);
}

// Revised Blossom Loop
int revised_blossom_loop(int ncount, int ecount, int *elist, double *x, int silent, CCrandstate *rstate) {
    int cutcount = 0, cut_added = 0;
    int outside = 0, num_loop = 10;
    CCtsp_lpcut_in *cuts = NULL;

    do {
        cut_added = 0;  // Reset cut_added for this outer loop iteration

        // Fast Blossoms
        printf("\nRunning Fast Blossoms...\n");
        CCtsp_fastblossom(&cuts, &cutcount, ncount, ecount, elist, x);
        if (cutcount > 0) {
            cut_added += cutcount;
            verify_and_print_comb(cuts, ncount, ecount, elist, x);
            free_cuts(cuts);
            cuts = NULL;
        } else {
            printf("Fast Blossoms found no cuts.\n");
        }

        // Groetschel-Holland Fast Blossoms
        printf("\nRunning Groetschel-Holland Fast Blossoms...\n");
        CCtsp_ghfastblossom(&cuts, &cutcount, ncount, ecount, elist, x);
        if (cutcount > 0) {
            cut_added += cutcount;
            verify_and_print_comb(cuts, ncount, ecount, elist, x);
            free_cuts(cuts);
            cuts = NULL;
        } else {
            printf("GH Fast Blossoms found no cuts.\n");
        }

        // Exact Blossoms
        printf("\nRunning Exact Blossoms...\n");
        CCtsp_exactblossom(&cuts, &cutcount, ncount, ecount, elist, x, rstate);
        if (cutcount > 0) {
            cut_added += cutcount;
            verify_and_print_comb(cuts, ncount, ecount, elist, x);
            free_cuts(cuts);
            cuts = NULL;
        } else {
            printf("Exact Blossoms found no cuts.\n");
        }

    } while (cut_added > 0 && ++outside < num_loop);  // Continue if cuts were added

    if (cuts) free_cuts(cuts);
    return 0;  // Return 0 explicitly since rval is no longer tracked
}


// Verify and Print Comb Details
void verify_and_print_comb(CCtsp_lpcut_in *cuts, int ncount, int ecount, int *elist, double *x) {
    CCtsp_lpcut_in *current_cut = cuts;
    int comb_index = 0;

    while (current_cut) {
        printf("\nInspecting Comb %d:\n", ++comb_index);
        double delta_H = 0.0, delta_T_sum = 0.0;
        int num_teeth = current_cut->cliquecount - 1;

        // Print the handle
        printf("Handle: ");
        CCtsp_lpclique *handle = &(current_cut->cliques[0]);
        for (int s = 0; s < handle->segcount; s++) {
            printf("[%d, %d] ", handle->nodes[s].lo, handle->nodes[s].hi);
        }
        printf("\n");
        // Print the teeth
        for (int t = 1; t <= num_teeth; t++) {
            printf("Tooth %d: ", t);
            CCtsp_lpclique *tooth = &(current_cut->cliques[t]);
            for (int s = 0; s < tooth->segcount; s++) {
                printf("[%d, %d] ", tooth->nodes[s].lo, tooth->nodes[s].hi);
            }
        }
        printf("\n");
        // Calculate delta_H and delta_T_sum
        for (int i = 0; i < ecount; i++) {
            int u = elist[2 * i];
            int v = elist[2 * i + 1];
            for (int t = 0; t < current_cut->cliquecount; t++) {
                CCtsp_lpclique *clique = &(current_cut->cliques[t]);
                int in_handle = 0, in_teeth = 0;
                for (int s = 0; s < clique->segcount; s++) {
                    if (clique->nodes[s].lo <= u && u <= clique->nodes[s].hi) in_handle++;
                    if (clique->nodes[s].lo <= v && v <= clique->nodes[s].hi) in_teeth++;
                }
                if ((in_handle ^ in_teeth)) {
                    if (t == 0) delta_H += x[i];
                    else delta_T_sum += x[i];
                }
            }
        }

        // Print the inequality
        printf("LHS = %.4f, RHS = %d\n", delta_H + delta_T_sum, 3 * num_teeth + 1);
        if (delta_H + delta_T_sum < 3 * num_teeth + 1) {
            printf("  Comb inequality is violated.\n");
        } else {
            printf("  Comb inequality is satisfied.\n");
            exit(1);
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
    printf("Freeing cuts completed\n");
}
