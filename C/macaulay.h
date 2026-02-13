#include <stdio.h>
#include <stdlib.h>
#include <flint/flint.h>
#include <flint/nmod_mpoly.h>
#include <m4ri/m4ri.h>

void macaulay_matrix(mzd_t *Mac, nmod_mpoly_t **Mx, unsigned int s, unsigned int limit, nmod_mpoly_ctx_t ctx);