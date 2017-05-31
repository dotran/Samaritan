#include "toolkit.h"
void wfg4(individual_real *ind)
{
    double *xreal, *obj;
    int size;
    int i;
    Degenerate = 0;

    WFG_ini();

    obj   = ind->obj;
    xreal = ind->xreal;
    size = number_variable;
    size = WFG4_t1(xreal,size,wfg_temp);
    size = WFG2_t3(wfg_temp,size,wfg_K,number_objective,wfg_temp);
    WFG4_shape(wfg_temp,size,obj);

    WFG_free();


}
void wfg42(individual_real *ind)
{
    double *xreal, *obj;
    int size;
    int i;
    Degenerate = 0;

    WFG_ini();

    obj   = ind->obj;
    xreal = ind->xreal;
    size = number_variable;
    size = WFG4_t1(xreal,size,wfg_temp);
    size = WFG2_t3(wfg_temp,size,wfg_K,number_objective,wfg_temp);
    WFG42_shape(wfg_temp,size,obj);

    WFG_free();


}


void wfg43(individual_real *ind)
{
    double *xreal, *obj;
    int size;
    int i;
    Degenerate = 0;

    WFG_ini();

    obj   = ind->obj;
    xreal = ind->xreal;
    size = number_variable;
    size = WFG4_t1(xreal,size,wfg_temp);
    size = WFG2_t3(wfg_temp,size,wfg_K,number_objective,wfg_temp);
    WFG43_shape(wfg_temp,size,obj);

    WFG_free();


}

void wfg44(individual_real *ind)
{
    double *xreal, *obj;
    int size;
    int i;
    Degenerate = 0;

    WFG_ini();

    obj   = ind->obj;
    xreal = ind->xreal;
    size = number_variable;
    size = WFG4_t1(xreal,size,wfg_temp);
    size = WFG2_t3(wfg_temp,size,wfg_K,number_objective,wfg_temp);
    WFG44_shape(wfg_temp,size,obj);

    WFG_free();


}

void wfg45(individual_real *ind)
{
    double *xreal, *obj;
    int size;
    int i;
    Degenerate = 0;

    WFG_ini();

    obj   = ind->obj;
    xreal = ind->xreal;
    size = number_variable;
    size = WFG4_t1(xreal,size,wfg_temp);
    size = WFG2_t3(wfg_temp,size,wfg_K,number_objective,wfg_temp);
    WFG45_shape(wfg_temp,size,obj);

    WFG_free();


}


void wfg46(individual_real *ind)
{
    double *xreal, *obj;
    int size;
    int i;
    Degenerate = 0;

    WFG_ini();

    obj   = ind->obj;
    xreal = ind->xreal;
    size = number_variable;
    size = WFG4_t1(xreal,size,wfg_temp);
    size = WFG2_t3(wfg_temp,size,wfg_K,number_objective,wfg_temp);
    WFG46_shape(wfg_temp,size,obj);

    WFG_free();


}


void wfg47(individual_real *ind)
{
    double *xreal, *obj;
    int size;
    int i;
    Degenerate = 0;

    WFG_ini();

    obj   = ind->obj;
    xreal = ind->xreal;
    size = number_variable;
    size = WFG4_t1(xreal,size,wfg_temp);
    size = WFG2_t3(wfg_temp,size,wfg_K,number_objective,wfg_temp);
    WFG47_shape(wfg_temp,size,obj);

    WFG_free();


}


void wfg48(individual_real *ind)
{
    double *xreal, *obj;
    int size;
    int i;
    Degenerate = 0;

    WFG_ini();

    obj   = ind->obj;
    xreal = ind->xreal;
    size = number_variable;
    size = WFG4_t1(xreal,size,wfg_temp);
    size = WFG2_t3(wfg_temp,size,wfg_K,number_objective,wfg_temp);
    WFG48_shape(wfg_temp,size,obj);

    WFG_free();


}