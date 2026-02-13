#include "pfaffian.h"

nmod_mpoly_t **generate_matrix(matrix_t* C_mat, int s, nmod_mpoly_ctx_t ctx)
{
    unsigned int nvars = ctx->minfo->nvars;
    unsigned int N = nvars + 1;
    nmod_mpoly_t **Mx = (nmod_mpoly_t**) malloc(s * sizeof(nmod_mpoly_t*));
    for (int i = 0; i < s; i++)
    {
        Mx[i] = (nmod_mpoly_t*) malloc(s * sizeof(nmod_mpoly_t));
        for (int j = 0; j < s; j++)
            nmod_mpoly_init(Mx[i][j], ctx);
    }

    // Generate locally the variables x_1,...,x_N
    nmod_mpoly_t vars[nvars];
    for (int i = 0; i < nvars; i++)
    {
        nmod_mpoly_init(vars[i], ctx);
        nmod_mpoly_gen(vars[i], i, ctx);
    }

    // compute x_1 M_1 + ... + x_{N-1} M_{N-1} + M_N

    nmod_mpoly_t tmp;
    nmod_mpoly_init(tmp, ctx);

    for (int i = 0; i < s; i++)
    {
        for (int j = 0; j < s; j++)
        {
            nmod_mpoly_zero(Mx[i][j], ctx);
            for (int k = 0; k < N-1; k++)
            {
                /* tmp = M_k[i,j] * x_k */
                nmod_mpoly_scalar_mul_ui(tmp, vars[k], C_mat[k][i][j], ctx);

                /* Mx[i,j] = Mx[i,j] + tmp */
                nmod_mpoly_add(Mx[i][j], Mx[i][j], tmp, ctx);
            }

            /* Add C_mat[N-1][i][j] */
            nmod_mpoly_add_ui(Mx[i][j], Mx[i][j], C_mat[N-1][i][j], ctx);
        }
    }



    nmod_mpoly_clear(tmp, ctx);
    for (int i = 0; i < nvars; i++)
        nmod_mpoly_clear(vars[i], ctx);
    
    return Mx;
}




nmod_mpoly_t *Pfaffians(nmod_mpoly_t **Mx, int s, nmod_mpoly_ctx_t ctx)
{
    int N_pf = (s * (s-1) * (s-2) * (s-3)) / 24;
    nmod_mpoly_t *pfs = (nmod_mpoly_t *) malloc(N_pf * sizeof(nmod_mpoly_t));

    nmod_mpoly_t term_1, term_2, term_3;
    nmod_mpoly_init(term_1, ctx);
    nmod_mpoly_init(term_2, ctx);
    nmod_mpoly_init(term_3, ctx);

    int count = 0;
    for (int i = 0; i < s; i++)
    {
        for (int j = i+1; j < s; j++)
        {
            for (int k = j+1; k < s; k++)
            {
                for (int l = k+1; l < s; l++)
                {
                    nmod_mpoly_init(pfs[count], ctx);

                    nmod_mpoly_mul(term_1, Mx[i][j], Mx[k][l], ctx);
                    nmod_mpoly_mul(term_2, Mx[i][k], Mx[j][l], ctx);
                    nmod_mpoly_mul(term_3, Mx[i][l], Mx[j][k], ctx);

                    nmod_mpoly_add(pfs[count], term_1, term_2, ctx);
                    nmod_mpoly_add(pfs[count], pfs[count], term_3, ctx);

                    count++;
                }
            }
        }
    }

    nmod_mpoly_clear(term_1, ctx);
    nmod_mpoly_clear(term_2, ctx);
    nmod_mpoly_clear(term_3, ctx);

    return pfs;
}