/* psl.c

Provide a library of gsl functions.

Features:
- request calculation via message
- function lookup

Author: shakfu
Repo: https://github.com/shakfu/pd-psl.git

*/
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_clausen.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_debye.h>
#include <math.h>

#include "m_pd.h"

/*
 * function lookup infratructure
 * ---------------------------------------------------------------------------
 */

// rudimentary str -> int conversion
// unsigned long hash(const char *str) {
//     unsigned int h = 0;
//     int c;

//     while ((c = *str++)) h += c;

//     return h;
// }

unsigned long hash(const char *str) {
    unsigned int h = 0;
    int c;

    while ((c = *str++)) 
        h += (h << 1) + c;

    return h;
}


// unsigned long hash(const char *str)
// {
//     unsigned long h = 5381;
//     int c;

//     while ((c = *str++))
//         h = ((h << 5) + h) + c; /* h * 33 + c */

//     return h;
// }



enum FUNC {
    LOG1P = 12931,
    EXPM1 = 12805,
    HYPOT = 13148,
    HYPOT3 = 39495,
    ACOSH = 11978,
    ASINH = 12341,
    ATANH = 12296,
    LDEXP = 12829,
    POW_INT = 122216,
    POW_2 = 13475,
    POW_3 = 13476,
    POW_4 = 13477,
    POW_5 = 13478,
    POW_6 = 13479,
    POW_7 = 13480,
    POW_8 = 13481,
    POW_9 = 13482,
    RANDO = 13254,
    AIRY_AI = 109980,
    AIRY_BI = 109983,
    BESSEL = 36546,
    CLAUSEN = 110879,
    DAWSON = 36848,
    DEBYE_1 = 109891,
    DEBYE_2 = 109892,
    DEBYE_3 = 109893,
    DEBYE_4 = 109894,
};

/*
 * forward declarations / prototypes
 * ---------------------------------------------------------------------------
 */

typedef struct _psl t_psl;

void select_default_function(t_psl *x, t_symbol *s);

/*
 * psl class objects
 * ---------------------------------------------------------------------------
 */

static t_class *psl_class;

/*
 * psl class struct (data-space)
 * ---------------------------------------------------------------------------
 */

typedef void (*unary_func)(t_psl *, t_floatarg);
typedef void (*binary_func)(t_psl *, t_floatarg, t_floatarg);
typedef void (*tri_func)(t_psl *, t_floatarg, t_floatarg, t_floatarg);

typedef struct _psl {
    t_object x_obj;

    // assigned function
    t_symbol *func_name;
    int nargs;

    // function slots
    unary_func ufunc;
    binary_func bfunc;
    tri_func tfunc;

    // private vars

    // outlets
    t_outlet *out_f;
} t_psl;

/*
 * psl class methods (operation-space)
 * ---------------------------------------------------------------------------
 */

// typed-methods
void psl_float(t_psl *x, t_floatarg f) {
    post("psl_float: %f", f);
    x->ufunc(x, f);
}

void psl_list(t_psl *x, t_symbol *s, int argc, t_atom *argv) {

    if (s == gensym("list")) {
        post("s: list");

        if (argc == 0) {
            return;
        }

        if (argc == 1) {
            if (argv->a_type == A_FLOAT) {
                float f = atom_getfloat(argv);
                post("got float: %f", f);
                x->ufunc(x, f);
                return;
            }
        }

        if (argc == 2) {
            char buf[1000];
            for (int i = 0; i < argc; i++) {
                atom_string((argv+i), buf, 1000);
                post("arg+%i: %s", i, buf);
            }
            if (argv->a_type == A_FLOAT && (argv + 1)->a_type == A_FLOAT) {
                float f1 = atom_getfloat(argv+0);
                float f2 = atom_getfloat(argv+1);
                post("f(%.2f, %.2f)", f1, f2);
                x->bfunc(x, f1, f2);
                return;
            }
        }

        if (argc == 3) {
            char buf[1000];
            for (int i = 0; i < argc; i++) {
                atom_string((argv+i), buf, 1000);
                post("arg+%i: %s", i, buf);
            }
            if (argv->a_type == A_FLOAT && (argv+1)->a_type == A_FLOAT && (argv+2)->a_type == A_FLOAT) {
                float f1 = atom_getfloat(argv+0);
                float f2 = atom_getfloat(argv+1);
                float f3 = atom_getfloat(argv+2);                
                post("f(%.2f, %.2f, %.2f)", f1, f2, f3);
                x->tfunc(x, f1, f2, f3);
                return;
            }
        }
    }

    post("list body");
    return;

    // error:
    //     pd_error(x, "psl_list error: incorrect arg type");
}

// message-methods
void psl_rando(t_psl *x, t_floatarg n, t_floatarg seed) {
    post("rando: n:%.2f seed:%.2f", n, seed);
    gsl_rng_env_setup();

    int argc = (int)n;
    int _seed = (int)seed;

    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, _seed);

    t_atom av[argc];
    for (int i = 0; i < argc; i++) {
        SETFLOAT(av + i, gsl_rng_uniform(r));
    }
    outlet_list(x->out_f, gensym("list"), argc, av);

    gsl_rng_free(r);
}

void psl_airy_ai(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_sf_airy_Ai(f, GSL_PREC_APPROX));
}

void psl_airy_bi(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_sf_airy_Bi(f, GSL_PREC_APPROX));
}

void psl_log1p(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_log1p(f));
}

void psl_expm1(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_expm1(f));
}

void psl_hypot(t_psl *x, t_floatarg f1, t_floatarg f2) {
    outlet_float(x->out_f, gsl_hypot(f1, f2));
}

void psl_hypot3(t_psl *x, t_floatarg f1, t_floatarg f2, t_floatarg f3) {
    outlet_float(x->out_f, gsl_hypot3(f1, f2, f3));
}

void psl_acosh(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_acosh(f));
}

void psl_asinh(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_asinh(f));
}

void psl_atanh(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_atanh(f));
}

void psl_ldexp(t_psl *x, t_floatarg f1, t_floatarg f2) {
    outlet_float(x->out_f, gsl_ldexp(f1, f2));
}

void psl_pow_int(t_psl *x, t_floatarg f1, t_floatarg f2) {
    outlet_float(x->out_f, gsl_pow_int(f1, f2));
}

void psl_pow_2(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_pow_2(f));
}

void psl_pow_3(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_pow_3(f));
}

void psl_pow_4(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_pow_4(f));
}

void psl_pow_5(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_pow_5(f));
}

void psl_pow_6(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_pow_6(f));
}

void psl_pow_7(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_pow_7(f));
}

void psl_pow_8(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_pow_8(f));
}

void psl_pow_9(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_pow_9(f));
}

void psl_bessel(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_sf_bessel_J0(f));
}

void psl_clausen(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_sf_clausen(f));
}

void psl_dawson(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_sf_dawson(f));
}

void psl_debye_1(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_sf_debye_1(f));
}

void psl_debye_2(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_sf_debye_2(f));
}

void psl_debye_3(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_sf_debye_3(f));
}

void psl_debye_4(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_sf_debye_4(f));
}


/*
 * function selection
 * ---------------------------------------------------------------------------
 */

// set default function from symbol
void select_default_function(t_psl *x, t_symbol *s) {
    x->func_name = s;
    post("func %s selected", s->s_name);

    switch (hash(s->s_name)) {
        case LOG1P:
            x->nargs = 1;
            x->ufunc = &psl_log1p;
            break;
        case EXPM1:
            x->nargs = 1;
            x->ufunc = &psl_expm1;
            break;
        case HYPOT:
            x->nargs = 2;
            x->bfunc = &psl_hypot;
            break;
        case HYPOT3:
            x->nargs = 3;
            x->tfunc = &psl_hypot3;
            break;
        case ACOSH:
            x->nargs = 1;
            x->ufunc = &psl_acosh;
            break;
        case ASINH:
            x->nargs = 1;
            x->ufunc = &psl_asinh;
            break;
        case ATANH:
            x->nargs = 1;
            x->ufunc = &psl_atanh;
            break;
        case LDEXP:
            x->nargs = 2;
            x->bfunc = &psl_ldexp;
            break;
        case POW_INT:
            x->nargs = 2;
            x->bfunc = &psl_pow_int;
            break;
        case POW_2:
            x->nargs = 1;
            x->ufunc = &psl_pow_2;
            break;
        case POW_3:
            x->nargs = 1;
            x->ufunc = &psl_pow_3;
            break;
        case POW_4:
            x->nargs = 1;
            x->ufunc = &psl_pow_4;
            break;
        case POW_5:
            x->nargs = 1;
            x->ufunc = &psl_pow_5;
            break;
        case POW_6:
            x->nargs = 1;
            x->ufunc = &psl_pow_6;
            break;
        case POW_7:
            x->nargs = 1;
            x->ufunc = &psl_pow_7;
            break;
        case POW_8:
            x->nargs = 1;
            x->ufunc = &psl_pow_8;
            break;
        case POW_9:
            x->nargs = 1;
            x->ufunc = &psl_pow_9;
            break;
        case RANDO:
            x->nargs = 2;
            x->bfunc = &psl_rando;
            break;
        case AIRY_AI:
            x->nargs = 1;
            x->ufunc = &psl_airy_ai;
            break;
        case AIRY_BI:
            x->nargs = 1;
            x->ufunc = &psl_airy_bi;
            break;
        case BESSEL:
            x->nargs = 1;
            x->ufunc = &psl_bessel;
            break;
        case CLAUSEN:
            x->nargs = 1;
            x->ufunc = &psl_clausen;
            break;
        case DAWSON:
            x->nargs = 1;
            x->ufunc = &psl_dawson;
            break;
        case DEBYE_1:
            x->nargs = 1;
            x->ufunc = &psl_debye_1;
            break;
        case DEBYE_2:
            x->nargs = 1;
            x->ufunc = &psl_debye_2;
            break;
        case DEBYE_3:
            x->nargs = 1;
            x->ufunc = &psl_debye_3;
            break;
        case DEBYE_4:
            x->nargs = 1;
            x->ufunc = &psl_debye_4;
            break;
        default:
            post("func selection failed, reverting to defaults");
            break;
   }
}

/*
 * psl class constructor
 * ---------------------------------------------------------------------------
 */

void *psl_new(t_symbol *s) {
    t_psl *x = (t_psl *)pd_new(psl_class);

    // initialize variables
    x->nargs = 1;
    x->ufunc = &psl_bessel;

    select_default_function(x, s);

    // create inlets

    // initialize outlets
    x->out_f = outlet_new(&x->x_obj, &s_float);

    return (void *)x;
}

/*
 * psl class setup
 * ---------------------------------------------------------------------------
 */

void psl_setup(void) {
    psl_class = class_new(gensym("psl"), (t_newmethod)psl_new,
                          0,  // destructor
                          sizeof(t_psl), CLASS_DEFAULT, A_DEFSYMBOL, 0);

    // typed methods
    class_addfloat(psl_class, psl_float);
    class_addlist(psl_class, psl_list);

    // message methods

    // functions
    class_addmethod(psl_class, (t_method)psl_log1p,  gensym("log1p"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_expm1,  gensym("expm1"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_hypot,  gensym("hypot"), A_DEFFLOAT, A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_hypot3,  gensym("hypot3"), A_DEFFLOAT, A_DEFFLOAT, A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_acosh,  gensym("acosh"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_asinh,  gensym("asinh"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_atanh,  gensym("atanh"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_ldexp,  gensym("ldexp"), A_DEFFLOAT, A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_pow_int,  gensym("pow_int"), A_DEFFLOAT, A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_pow_2,  gensym("pow_2"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_pow_3,  gensym("pow_3"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_pow_4,  gensym("pow_4"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_pow_5,  gensym("pow_5"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_pow_6,  gensym("pow_6"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_pow_7,  gensym("pow_7"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_pow_8,  gensym("pow_8"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_pow_9,  gensym("pow_9"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_rando,  gensym("rando"), A_DEFFLOAT, A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_airy_ai,  gensym("airy_ai"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_airy_bi,  gensym("airy_bi"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_bessel,  gensym("bessel"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_clausen,  gensym("clausen"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_dawson,  gensym("dawson"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_debye_1,  gensym("debye_1"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_debye_2,  gensym("debye_2"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_debye_3,  gensym("debye_3"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_debye_4,  gensym("debye_4"), A_DEFFLOAT, 0);

    // set name of default help file
    class_sethelpsymbol(psl_class, gensym("help-psl"));
}
