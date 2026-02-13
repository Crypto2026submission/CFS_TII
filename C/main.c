#include <stdio.h>
#include <flint/flint.h>
#include <flint/nmod_mpoly.h>
#include <time.h>
#include <m4ri/m4ri.h>

#include "pfaffian.h"
#include "macaulay.h"
#include "read_cmat.h"

int main(void)
{
    m4ri_init();

    unsigned int nvars, s;
    matrix_t *Cmat = read_file("Cmat.txt", &s, &nvars);

    nvars--;
    printf("nvars = %d\n", nvars);
    /* Polynomial context over GF(2) */
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_ctx_init(ctx, nvars, ORD_DEGLEX, 2);

    nmod_mpoly_t **Mx = generate_matrix(Cmat, s, ctx);

    unsigned int N_pfs = ((s) * (s - 1) * (s - 2) * (s - 3)) / 24;
    unsigned int N_mons = ((nvars + 1) * (nvars + 2)) / 2;

    /* limit encodes the number of pfaffians that are computed. Change it if need be */
    unsigned int limit = 105000;
    unsigned int N_rows;

    if (limit < N_pfs)
        N_rows = limit;
    
    else 
        N_rows = N_pfs;

    printf("Computing %d pfaffians and the Macaulay matrix...\n", N_rows);
    
    clock_t T0 = clock();

    mzd_t *Mac = mzd_init(N_rows, N_mons);
    macaulay_matrix(Mac, Mx, s, limit, ctx);

    clock_t T1 = clock();
    printf("Time: %f s\n", (double)(T1 - T0) / CLOCKS_PER_SEC);

    puts("Computing reduced row echelon form...");
    slong rank = mzd_echelonize(Mac, 0);
    clock_t T2 = clock();
    printf("Time: %f s\n", (double)(T2 - T1) / CLOCKS_PER_SEC);

    printf("Macaulay matrix has rank %ld\n", rank);

    /* Print the expected linear equations */
    nmod_mpoly_t linear;
    nmod_mpoly_init(linear, ctx);

    nmod_mpoly_t vars[nvars];
    for (int i = 0; i < nvars; i++)
    {
        nmod_mpoly_init(vars[i], ctx);
        nmod_mpoly_gen(vars[i], i, ctx);
    }

    unsigned int row = rank-1;
    unsigned int Nmons = nvars * (nvars + 1) / 2;

    while (1)
    {
        int check = 1;
        for (int i = 0; i < Nmons; i++)
        {
            check = check && !mzd_read_bit(Mac,row,i);
        }
        
        if (!check)
            break;

        nmod_mpoly_zero(linear, ctx);

        for (int i = 0; i < nvars; i++)
        {
            if (mzd_read_bit(Mac,row,Nmons+i))
            {
                nmod_mpoly_add(linear, linear, vars[i], ctx);
            }
        }
        nmod_mpoly_add_ui(linear, linear, mzd_read_bit(Mac,row,N_mons-1), ctx);

        nmod_mpoly_print_pretty(linear, NULL, ctx);
        puts("");
        row--;
    }

    printf("%d linear equations computed\n", rank - row - 1);

    /* Cleanup */
    printf("Total time: %f s\n", (double)(T2 - T0) / CLOCKS_PER_SEC);
    for (int i = 0; i < s; ++i)
    {
        for (int j = 0; j < s; ++j)
        {
            nmod_mpoly_clear(Mx[i][j], ctx);
        }
        free(Mx[i]);
    }
    free(Mx);

    nmod_mpoly_clear(linear, ctx);

    for (int i = 0; i < nvars; i++)
    {
        nmod_mpoly_clear(vars[i], ctx);
    }

    nmod_mpoly_ctx_clear(ctx);

    for (int k = 0; k < nvars+1; k++)
    {
        for (int i = 0; i < s; i++)
        {
            free(Cmat[k][i]);
        }
        free(Cmat[k]);
    }
    free(Cmat);

    mzd_free(Mac);
    m4ri_fini();

    return 0;
}
