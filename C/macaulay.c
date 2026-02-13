#include "macaulay.h"

void macaulay_matrix(mzd_t *Mac, nmod_mpoly_t **Mx, unsigned int s, unsigned int limit, nmod_mpoly_ctx_t ctx)
{
    unsigned int nvars  = ctx->minfo->nvars;
    unsigned int N_mons = ((nvars + 1)*(nvars + 2))/2;
    unsigned int N_eqs  = (s * (s-1) * (s-2) * (s-3)) / 24;
    unsigned int N_rows;

    if (limit < N_eqs)
    {
        N_rows = limit;
    }

    else
    {
        N_rows = N_eqs;
    }

    if (Mac->nrows != N_rows || Mac->ncols != N_mons)
        exit(EXIT_FAILURE);

    /* Create a second context in odd characteristics */
    nmod_mpoly_ctx_t ctx2;
    nmod_mpoly_ctx_init(ctx2, nvars, ORD_DEGLEX, 3);

    /* Create a vector containing all variables */
    nmod_mpoly_t vars[nvars];
    for (int i = 0; i < nvars; ++i)
    {
        nmod_mpoly_init(vars[i], ctx2);
        nmod_mpoly_gen(vars[i], i, ctx2);
    }

    /* Create a vector containing all monomials of degree at most 2 */
    puts("Creating sum");
    nmod_mpoly_t sum;
    nmod_mpoly_init2(sum, N_mons, ctx2);

    // Add 1
    nmod_mpoly_add_ui(sum, sum, 1, ctx2);
    puts("sum = 1");

    // Add all variables
    for (int i = 0; i < nvars; ++i)
    {
        nmod_mpoly_add(sum, sum, vars[i], ctx2);
    }
    puts("Variables added");

    nmod_mpoly_mul(sum, sum, sum,ctx2);


    puts("All monomials added");


    nmod_mpoly_t pf;
    nmod_mpoly_init2(pf, N_mons, ctx);

    nmod_mpoly_t term_1, term_2, term_3;
    nmod_mpoly_init2(term_1, N_mons, ctx);
    nmod_mpoly_init2(term_2, N_mons, ctx);
    nmod_mpoly_init2(term_3, N_mons, ctx);

    int count;
    int row = 0;
    int N = mpoly_words_per_exp(sum->bits, ctx->minfo);
    for (int i = 0; i < s; i++)
    {
        for (int j = i+1; j < s; j++)
        {
            for (int k = j+1; k < s; k++)
            {
                for (int l = k+1; l < s; l++)
                {

                    nmod_mpoly_zero(pf, ctx);
                    // printf("i,j,k,l = %d,%d,%d,%d\n", i,j,k,l);

                    nmod_mpoly_mul(term_1, Mx[i][j], Mx[k][l], ctx);
                    nmod_mpoly_mul(term_2, Mx[i][k], Mx[j][l], ctx);
                    nmod_mpoly_mul(term_3, Mx[i][l], Mx[j][k], ctx);

                    nmod_mpoly_add(pf, term_1, term_2, ctx);
                    nmod_mpoly_add(pf, pf, term_3, ctx);

                    // Now that pf is computed, write the corresponding row in Mac
                    count = 0;

                    for (int m = 0; m < pf->length; m++)
                    {
                        while (!(mpoly_monomial_equal(pf->exps + N*m, sum->exps + N*count, N)))
                            count++;
                        mzd_write_bit(Mac, row, count, 1);
                    }
                    row++;

                    if (row % 10000 == 0)
                    {
                        printf("%f %%\n", (((double) 100 * row) / N_rows));
                    }

                    if (row == limit)
                        goto cleanup;
                }
            }
        }
    }

    cleanup:
        nmod_mpoly_clear(term_1, ctx);
        nmod_mpoly_clear(term_2, ctx);
        nmod_mpoly_clear(term_3, ctx);

        nmod_mpoly_clear(pf, ctx);

        for (int i = 0; i < nvars; ++i)
            nmod_mpoly_clear(vars[i], ctx2);
        
        nmod_mpoly_clear(sum, ctx2);
        nmod_mpoly_ctx_clear(ctx2);
}