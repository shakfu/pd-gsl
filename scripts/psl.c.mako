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
    % for f in funcs:
    ${f.name.upper()} = ${f.hashed},
    % endfor
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

% for f in funcs:
% if f.name not in skip:
% if f.nargs == 1:
void psl_${f.name}(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_${f.fullname}(f));
}

% endif
% if f.nargs == 2:
void psl_${f.name}(t_psl *x, t_floatarg f1, t_floatarg f2) {
    outlet_float(x->out_f, gsl_${f.fullname}(f1, f2));
}

% endif
% if f.nargs == 3:
void psl_${f.name}(t_psl *x, t_floatarg f1, t_floatarg f2, t_floatarg f3) {
    outlet_float(x->out_f, gsl_${f.fullname}(f1, f2, f3));
}

% endif
% endif
% endfor

/*
 * function selection
 * ---------------------------------------------------------------------------
 */

// set default function from symbol
void select_default_function(t_psl *x, t_symbol *s) {
    x->func_name = s;
    post("func %s selected", s->s_name);

    switch (hash(s->s_name)) {
        % for f in funcs:
        case ${f.name.upper()}:
            x->nargs = ${f.nargs};
            x->${f.ftype} = &psl_${f.name};
            break;
        % endfor
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
    x->ufunc = &psl_bessel_j0;

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

    // class-addmethods
    % for f in funcs:
    class_addmethod(psl_class, (t_method)psl_${f.name},  gensym("${f.name}"), ${f.slots}, 0);
    % endfor

    // create alias
    class_addcreator((t_newmethod)psl_new, gensym("gsl"), A_DEFSYMBOL, 0);

    // set name of default help file
    class_sethelpsymbol(psl_class, gensym("help-psl"));
}
