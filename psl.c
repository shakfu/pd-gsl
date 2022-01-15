/* psl.c
////
Provide a library of gsl functions.

Features:
- request calculation via message
- function lookup
- inlets-on-demand

Author: shakfu
Repo: https://github.com/shakfu/pd-psl.git

*/
#include <math.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_clausen.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_debye.h>

#include "m_pd.h"
#include "tinyexpr.h"


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
    ADD = 1273,
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
    FCMP = 4084,
    AIRY_AI = 109980,
    AIRY_BI = 109983,
    BESSEL_J0 = 987963,
    BESSEL_J1 = 987964,
    BESSEL_JN = 988025,
    BESSEL_Y0 = 988008,
    BESSEL_Y1 = 988009,
    BESSEL_YN = 988070,
    BESSEL_I0 = 987960,
    BESSEL_I1 = 987961,
    BESSEL_IN = 988022,
    CLAUSEN = 110879,
    DAWSON = 36848,
    DEBYE_1 = 109891,
    DEBYE_2 = 109892,
    DEBYE_3 = 109893,
    DEBYE_4 = 109894,
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

    // for expression
    char expr_buffer[MAXPDSTRING];

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

    // atom_post("psl_list: ", argc, argv);

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



void psl_symbol(t_psl *x, t_symbol *s) {
    post("s: %s", s->s_name);

    // local buffer
    int length = strlen(s->s_name);
    char *buf = (char *)malloc(length * sizeof(char));
    strcpy(buf, s->s_name);
    post("buf: %s", buf);

    // clear expr_buffer
    memset(x->expr_buffer, 0, MAXPDSTRING);

    int j = 0;
    for (int i = 0; i < length; i++) {
        // remove escape `\` required for commas
        if (buf[i] != '\\') {
            x->expr_buffer[j++] = buf[i];
        } else if (x->expr_buffer[j - 1] == ' ') {
            j--;
        }
    }
    x->expr_buffer[length] = '\0';
    free(buf);

    post("x->expr_buffer: %s", x->expr_buffer);

    te_variable vars[] = {
        {"hypot", gsl_hypot, TE_FUNCTION2, NULL} /* TE_FUNCTION2 used because my_sum takes two arguments. */
    };

    te_expr *expr = te_compile(x->expr_buffer, vars, 2, 0);
    const double res = te_eval(expr);
    te_free(expr);
    outlet_float(x->out_f, res);
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

void psl_fcmp(t_psl *x, t_floatarg f1, t_floatarg f2, t_floatarg f3) {
    outlet_float(x->out_f, gsl_fcmp(f1, f2, f3));
}

void psl_bessel_j0(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_sf_bessel_J0(f));
}

void psl_bessel_j1(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_sf_bessel_J1(f));
}

void psl_bessel_jn(t_psl *x, t_floatarg f1, t_floatarg f2) {
    outlet_float(x->out_f, gsl_sf_bessel_Jn(f1, f2));
}

void psl_bessel_y0(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_sf_bessel_Y0(f));
}

void psl_bessel_y1(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_sf_bessel_Y1(f));
}

void psl_bessel_yn(t_psl *x, t_floatarg f1, t_floatarg f2) {
    outlet_float(x->out_f, gsl_sf_bessel_Yn(f1, f2));
}

void psl_bessel_i0(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_sf_bessel_I0(f));
}

void psl_bessel_i1(t_psl *x, t_floatarg f) {
    outlet_float(x->out_f, gsl_sf_bessel_I1(f));
}

void psl_bessel_in(t_psl *x, t_floatarg f1, t_floatarg f2) {
    outlet_float(x->out_f, gsl_sf_bessel_In(f1, f2));
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



// function selection
//---------------------------------------------------------------------------


// set default function from symbol
void select_default_function(t_psl *x, t_symbol *s) {
    x->func_name = s;
    post("func %s selected", s->s_name);

    switch (hash(s->s_name)) {
        case ADD:
            x->nargs = 2;
            x->bfunc = &psl_add;
            break;
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
        case FCMP:
            x->nargs = 3;
            x->tfunc = &psl_fcmp;
            break;
        case AIRY_AI:
            x->nargs = 1;
            x->ufunc = &psl_airy_ai;
            break;
        case AIRY_BI:
            x->nargs = 1;
            x->ufunc = &psl_airy_bi;
            break;
        case BESSEL_J0:
            x->nargs = 1;
            x->ufunc = &psl_bessel_j0;
            break;
        case BESSEL_J1:
            x->nargs = 1;
            x->ufunc = &psl_bessel_j1;
            break;
        case BESSEL_JN:
            x->nargs = 2;
            x->bfunc = &psl_bessel_jn;
            break;
        case BESSEL_Y0:
            x->nargs = 1;
            x->ufunc = &psl_bessel_y0;
            break;
        case BESSEL_Y1:
            x->nargs = 1;
            x->ufunc = &psl_bessel_y1;
            break;
        case BESSEL_YN:
            x->nargs = 2;
            x->bfunc = &psl_bessel_yn;
            break;
        case BESSEL_I0:
            x->nargs = 1;
            x->ufunc = &psl_bessel_i0;
            break;
        case BESSEL_I1:
            x->nargs = 1;
            x->ufunc = &psl_bessel_i1;
            break;
        case BESSEL_IN:
            x->nargs = 2;
            x->bfunc = &psl_bessel_in;
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
    class_addsymbol(psl_class, psl_symbol);


    // message methods

    // class-addmethods
    class_addmethod(psl_class, (t_method)psl_add,  gensym("add"), A_DEFFLOAT, A_DEFFLOAT, 0);
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
    class_addmethod(psl_class, (t_method)psl_fcmp,  gensym("fcmp"), A_DEFFLOAT, A_DEFFLOAT, A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_airy_ai,  gensym("airy_ai"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_airy_bi,  gensym("airy_bi"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_bessel_j0,  gensym("bessel_j0"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_bessel_j1,  gensym("bessel_j1"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_bessel_jn,  gensym("bessel_jn"), A_DEFFLOAT, A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_bessel_y0,  gensym("bessel_y0"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_bessel_y1,  gensym("bessel_y1"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_bessel_yn,  gensym("bessel_yn"), A_DEFFLOAT, A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_bessel_i0,  gensym("bessel_i0"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_bessel_i1,  gensym("bessel_i1"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_bessel_in,  gensym("bessel_in"), A_DEFFLOAT, A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_clausen,  gensym("clausen"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_dawson,  gensym("dawson"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_debye_1,  gensym("debye_1"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_debye_2,  gensym("debye_2"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_debye_3,  gensym("debye_3"), A_DEFFLOAT, 0);
    class_addmethod(psl_class, (t_method)psl_debye_4,  gensym("debye_4"), A_DEFFLOAT, 0);

    // create alias
    class_addcreator((t_newmethod)psl_new, gensym("gsl"), A_DEFSYMBOL, 0);

    // set name of default help file
    class_sethelpsymbol(psl_class, gensym("help-psl"));
}
