#include <gsl/gsl_rng.h>
#include <stdio.h>

/*
gsl_rng_mt19937
gsl_rng_ranlxs0
gsl_rng_ranlxs1
gsl_rng_ranlxs2
*/

void print_algos(void) {

    const gsl_rng_type **t, **t0;

    t0 = gsl_rng_types_setup();

    printf("Available generators:\n");

    for (t = t0; *t != 0; t++) {
        printf("%s\n", (*t)->name);
    }

}

int main(void) {
    int seed = 102;
    int n = 10;

    gsl_rng_env_setup();

    gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, seed);

    for (int i = 0; i < n; i++) {
        double u = gsl_rng_uniform(r);
        printf("%.5f\n", u);
    }

    gsl_rng_free(r);

    return 0;
}
