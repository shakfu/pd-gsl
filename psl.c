/* psl.c

An external which counts via a variable step and optionally between two limits.


Features:
- configurable integer counting
- can count in steps
- can count between a lower and upper bound

Author: gpt3
Repo: https://github.com/gpt3/psl.git

*/
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_clausen.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_debye.h>


#include "m_pd.h"


/*
 * utility functions
 * ---------------------------------------------------------------------------
 */

unsigned long hash(char *str)
{
    unsigned int hash = 0;
    int c;

    while (c = *str++)
        hash += c;

    return hash;
}

enum FUNC {
    BESSEL = 638,
    DAWSON = 652
};


/*
 * psl class object
 * ---------------------------------------------------------------------------
 */

static t_class *psl_class;


/*
 * psl class struct (data-space)
 * ---------------------------------------------------------------------------
 */

typedef struct _psl t_psl;

typedef void (*simplefunc)(t_psl *, t_floatarg);

typedef struct _psl {
    t_object x_obj;

    // parameters
    t_float result;

    // function slot
    simplefunc func;

    // private vars

    // outlets
    t_outlet *out_f;
} t_psl;





/*
 * psl class methods (operation-space)
 * ---------------------------------------------------------------------------
 */

// typed-methods
void psl_bang(t_psl *x)
{
    post("bang -> outpuet last result");
}


void psl_float(t_psl *x, t_floatarg f)
{
    x->func(x, f);
}

void psl_list(t_psl *x, t_symbol *s, int argc, t_atom *argv)
{
    post("list body");
}


// message-methods
void psl_rando(t_psl *x, t_floatarg n, t_floatarg seed)
{
    gsl_rng_env_setup();

    gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, seed);

    for (int i = 0; i < n; i++) {
        outlet_float(x->out_f, gsl_rng_uniform(r));
    }

    gsl_rng_free(r);    
}

void psl_airy_ai(t_psl *x, t_floatarg f)
{
    outlet_float(x->out_f, gsl_sf_airy_Ai(f, GSL_PREC_APPROX));
}

void psl_airy_bi(t_psl *x, t_floatarg f)
{
    outlet_float(x->out_f, gsl_sf_airy_Bi(f, GSL_PREC_APPROX));
}

void psl_bessel(t_psl *x, t_floatarg f)
{
    outlet_float(x->out_f, gsl_sf_bessel_J0(f));
}

void psl_clausen(t_psl *x, t_floatarg f)
{
    outlet_float(x->out_f, gsl_sf_clausen(f));
}

// Oone-sided Fourier–Laplace sine transform of the Gaussian function.
void psl_dawson(t_psl *x, t_floatarg f)
{
    outlet_float(x->out_f, gsl_sf_dawson(f));
}

void psl_debye_1(t_psl *x, t_floatarg f)
{
    outlet_float(x->out_f, gsl_sf_debye_1(f));
}

void psl_debye_2(t_psl *x, t_floatarg f)
{
    outlet_float(x->out_f, gsl_sf_debye_2(f));
}

void psl_debye_3(t_psl *x, t_floatarg f)
{
    outlet_float(x->out_f, gsl_sf_debye_3(f));
}

void psl_debye_4(t_psl *x, t_floatarg f)
{
    outlet_float(x->out_f, gsl_sf_debye_4(f));
}




/*
 * psl class constructor
 * ---------------------------------------------------------------------------
 */

void *psl_new(t_symbol *s, t_floatarg f)
{
    t_psl *x = (t_psl *)pd_new(psl_class);

    // initialize variables
    x->result = 0.0;

    x->func = &psl_bessel;

    switch (hash(s->s_name)) {
        case BESSEL:
            x->func = &psl_bessel;
            break;
        case DAWSON:
            x->func = &psl_dawson;
            break;
        default:
            break;
    }

    // create inlets
    // inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("list"), gensym("bound"));
    // floatinlet_new(&x->x_obj, &x->step);

    // initialize outlets
    x->out_f = outlet_new(&x->x_obj, &s_float);

    return (void *)x;
}


/*
 * psl class setup
 * ---------------------------------------------------------------------------
 */

void psl_setup(void) 
{
    psl_class = class_new(gensym("psl"),
                            (t_newmethod)psl_new,
                            0, // destructor
                            sizeof(t_psl),
                            CLASS_DEFAULT,
                            A_DEFSYMBOL,
                            A_DEFFLOAT,
                            0);

    // typed methods
    class_addbang(psl_class, psl_bang);
    class_addfloat(psl_class, psl_float);
    class_addlist(psl_class, psl_list);

    // message methods
    class_addmethod(psl_class, (t_method)psl_rando,   gensym("rando"),   A_DEFFLOAT, A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_bessel,  gensym("bessel"),  A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_airy_ai, gensym("airy_ai"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_airy_bi, gensym("airy_bi"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_clausen, gensym("clausen"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_dawson,  gensym("dawson"),  A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_debye_1, gensym("debye_1"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_debye_2, gensym("debye_2"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_debye_3, gensym("debye_3"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_debye_4, gensym("debye_4"), A_DEFFLOAT, 0);


    // set name of default help file
    class_sethelpsymbol(psl_class, gensym("help-psl"));
}

