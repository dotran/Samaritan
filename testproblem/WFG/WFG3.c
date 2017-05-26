#include "toolkit.h"
void wfg3(individual_real *ind)
{
    double *xreal, *obj;
    int size;
    int i;
    Degenerate = 1;

    WFG_ini();

    obj   = ind->obj;
    xreal = ind->xreal;
    size = number_variable;
    size = WFG1_t1(xreal,size,wfg_K,wfg_temp);
    size = WFG2_t2(wfg_temp,size,wfg_K,wfg_temp);
    size = WFG2_t3(wfg_temp,size,wfg_K,number_objective,wfg_temp);
    WFG3_shape(wfg_temp,size,obj);

    WFG_free();


}