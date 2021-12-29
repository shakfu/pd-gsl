/* psl.c

Provide a library of gsl functions.

Features:
- request calculation via message
- function lookup
- inlets-on-demand

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


// macros and defines
//  ---------------------------------------------------------------------------


#define MAX_ARGS 6
#define STR_BUF_SIZE 1000


// function lookup infratructure
// ---------------------------------------------------------------------------

unsigned long hash(const char *str) {
    unsigned int h = 0;
    int c;

    while ((c = *str++)) 
        h += (h << 1) + c;

    return h;
}


enum FUNC {
    % for f in funcs:
    ${f.name.upper()} = ${f.hashed},
    % endfor
};


// forward declarations / prototypes
// ---------------------------------------------------------------------------


typedef struct _psl t_psl;

void select_default_function(t_psl *x, t_symbol *s);


// psl class objects
// ---------------------------------------------------------------------------


static t_class *psl_class;

static t_class *psl_inlet_class;


// psl class struct (data-space)
// ---------------------------------------------------------------------------


typedef void (*unary_func)(t_psl *, t_floatarg);
typedef void (*binary_func)(t_psl *, t_floatarg, t_floatarg);
typedef void (*tri_func)(t_psl *, t_floatarg, t_floatarg, t_floatarg);


typedef struct _psl_inlet
{
    t_class *x_pd;  // minimal pd object.
    t_psl   *owner; // the owning object to forward inlet messages to.
    int     id;     // the number of this inlet.
} t_psl_inlet;


typedef struct _psl {
    t_object x_obj;

    // assigned function
    t_symbol *func_name;
    int nargs;

    // function slots
    unary_func ufunc;
    binary_func bfunc;
    tri_func tfunc;

    // param_array
    t_float *arg_array;

    // private vars

    // inlets
    int inlets;          // # of extra inlets in addition to default
    t_psl_inlet *ins;    // the inlets themselves

    // outlets
    t_outlet *out_f;
} t_psl;


// psl class methods (operation-space)
// ---------------------------------------------------------------------------


// typed-methods

void psl_bang(t_psl *x) {
    if (x->nargs == 1 && x->inlets == 0) {
        x->ufunc(x, x->arg_array[0]);
    }

    if (x->nargs == 2 && x->inlets == 1) {
        x->bfunc(x, x->arg_array[0], x->arg_array[1]);
    }

    if (x->nargs == 3 && x->inlets == 2) {
        x->tfunc(x, x->arg_array[0], x->arg_array[1], x->arg_array[2]);
    }

}

void psl_float(t_psl *x, t_floatarg f) {
    post("psl_float: %f", f);
    if (x->nargs > 0) {
        x->arg_array[0] = f;
        psl_bang(x);
    } else {
        post("nothing to do: no function selected.");
        outlet_float(x->out_f, f);
    }
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
            char buf[STR_BUF_SIZE];
            for (int i = 0; i < argc; i++) {
                atom_string((argv+i), buf, STR_BUF_SIZE);
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
            char buf[STR_BUF_SIZE];
            for (int i = 0; i < argc; i++) {
                atom_string((argv+i), buf, STR_BUF_SIZE);
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

void psl_add(t_psl *x, t_floatarg f1, t_floatarg f2) {
    outlet_float(x->out_f, f1+f2);
}

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


// function selection
//---------------------------------------------------------------------------


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


// psl-inlet funcs
// ---------------------------------------------------------------------------


static void psl_inlet_float(t_psl_inlet *x, float f)
{
    x->owner->arg_array[x->id+1] = f;
    // outlet_float(x->owner->out_f, x->id + f);
    post("x->owner->arg_array[x->id]: %.02f", x->owner->arg_array[x->id+1]);
    psl_bang(x->owner);
}


// psl class constructor
// ---------------------------------------------------------------------------


void *psl_new(t_symbol *s) {
    t_psl *x = (t_psl *)pd_new(psl_class);

    // initialize variables
    x->nargs = 0;
    x->inlets = 0;
    x->ufunc = NULL;
    x->bfunc = NULL;
    x->tfunc = NULL;

    select_default_function(x, s);
    // sets x->nargs to correct number

    // create inlets
    x->inlets = x->nargs - 1;
    x->ins = (t_psl_inlet *)getbytes(x->inlets * sizeof(*x->ins));
    x->arg_array = malloc(x->inlets * sizeof(float));
    
    for (int i=0; i < x->inlets; i++) {
        x->ins[i].x_pd = psl_inlet_class;
        x->ins[i].owner = x;
        x->ins[i].id = i;
        x->arg_array[i] = 0.0;
        inlet_new((t_object *)x, &(x->ins[i].x_pd), 0, 0);
    }


    // initialize outlets
    x->out_f = outlet_new(&x->x_obj, &s_float);

    return (void *)x;
}


// psl class destructor
// ---------------------------------------------------------------------------


// TODO: not sure if this is correct!
void psl_free(t_psl *x) {
    free(x->arg_array);
    freebytes(x->ins, x->inlets * sizeof(*x->ins));
    post("DONE");
}


// psl class setup
// ---------------------------------------------------------------------------


void psl_setup(void) {

    psl_inlet_class = class_new(gensym("psl-inlet"), 
                                0, 0, 
                                sizeof(t_psl_inlet),
                                CLASS_PD,
                                0);

    if (psl_inlet_class) {
        class_addfloat(psl_inlet_class, (t_method)psl_inlet_float);
    }

    psl_class = class_new(gensym("psl"),
                        (t_newmethod)psl_new,
                        (t_method)psl_free,  // destructor
                        sizeof(t_psl), 
                        CLASS_DEFAULT, 
                        A_DEFSYMBOL, 
                        0);

    // typed methods
    class_addbang(psl_class, psl_bang);
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
