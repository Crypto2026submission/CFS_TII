#include "read_cmat.h"

u_int8_t ***read_file(char *filename, unsigned int *s, unsigned int *N)
{
    FILE *file;
    if ((file = fopen(filename, "r")) == NULL)
    {
        fprintf(stderr, "Cannot open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    if (fscanf(file, "%d", s) != 1) {
        fprintf(stderr, "Error reading first number\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    if (fscanf(file, "%d", N) != 1) {
        fprintf(stderr, "Error reading second number\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    u_int8_t ***Cmat = (u_int8_t ***) malloc(*N * (sizeof(u_int8_t**)));
    for (int k = 0; k < *N; k++)
    {
        Cmat[k] = (u_int8_t **) malloc(*s * sizeof(u_int8_t*));
        for (int i = 0; i < *s; i++)
        {
            Cmat[k][i] = (u_int8_t *) malloc(*s * sizeof(u_int8_t));
        }
    }

    getc(file);

    // read the matrices
    for (int k = 0; k < *N; k++)
    {
        for (int i = 0; i < *s; i++)
        {
            u_int8_t row[*s+2];
            if (fgets((char *)row, *s+2, file) == NULL)
            {
                fprintf(stderr, "Problem while reading a line\n");
                fclose(file);
                exit(EXIT_FAILURE);
            }

            for (int j = 0; j < *s; j++)
            {
                Cmat[k][i][j] = row[j];
            }
        }
    }
    fclose(file);

    return Cmat;
}