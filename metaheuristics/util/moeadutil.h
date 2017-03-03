/*
 * moeadutil.h:
 *  This is the header file for the utility functions for moead.
 *
 * Authors:
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *  Ke Li <k.li@exeter.ac.uk>
 *
 * Institution:
 *  Computational Optimization and Data Analytics (CODA) Group @ University of Exeter
 *
 * Copyright (c) 2017 Renzhi Chen, Ke Li
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef SAMARITAN_MOEADUTIL_H
#define SAMARITAN_MOEADUTIL_H
#include <stdio.h>
#include "../../header/global.h"
#include "../../header/utility.h"
#include "../../header/print.h"
#include "../../header/rand.h"


void setweight( double* v, double unit, double sum, int objdim, int dim);
void initialize_uniform_weight ();
void initialize_neighborhood();
void moead_free();
void initialize_idealpoint(void * pop);
void update_ideal_point(individual_real * p);
void update_neighborhood(population_real* pop,individual_real *individual, int subProblemId, int neighborType);
double fitnessFunction(individual_real* individual, double* lambda);
#endif