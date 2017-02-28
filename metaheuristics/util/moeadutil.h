//
// Created by rxc332 on 17-2-27.
//
#ifndef SAMARITAN_MOEADUTIL_H
#define SAMARITAN_MOEADUTIL_H
#include <stdio.h>
#include "../../header/global.h"
#include "../../header/utility.h"
#include "../../header/print.h"
#include "../../header/rand.h"
typedef struct double_s
{
    int idx;
    double x;
}mysort;

void setweight( double* v, double unit, double sum, int objdim, int dim);
void initialize_uniform_weight ();
void initialize_neighborhood();
int sort_double_cmp (const void * a, const void * b);
void moead_free();
void initialize_idealpoint(void * pop);
void update_ideal_point(individual_real * p);
void random_permutation(int* perm, int size);
void update_neighborhood(population_real* pop,individual_real *individual, int subProblemId, int neighborType);
double fitnessFunction(individual_real* individual, double* lambda);
#endif