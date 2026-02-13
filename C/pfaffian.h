#include <flint/flint.h>
#include <flint/nmod_mpoly.h>

#include "define.h"

nmod_mpoly_t **generate_matrix(matrix_t* C_mat, int s, nmod_mpoly_ctx_t ctx);

nmod_mpoly_t *Pfaffians(nmod_mpoly_t **Mx, int s, nmod_mpoly_ctx_t ctx);