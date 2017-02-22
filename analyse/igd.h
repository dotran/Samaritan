//
// Created by rxc332 on 17-2-21.
//

#ifndef SAMARITAN_IGD_H
#define SAMARITAN_IGD_H
//
// Created by rxc332 on 17-2-20.
//
#include "analyse.h"
#include <float.h>
#include "../header/print.h"
void igd(void *ptr,int id);
void print_igd(char * file_name);
void print_global_igd(char *file_name);
double point_distance(double *a, double *b, int dim);

static double *record = NULL;
static double *record_all_run = NULL;

#endif //SAMARITAN_IGD_H
