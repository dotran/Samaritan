/*
 * global.h:
 *  This is the header file for global variables and data structures.
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

# ifndef Samaritan_GLOBAL_H
# define Samaritan_GLOBAL_H

# include <stdio.h>
# include <stdlib.h>
# include <stdarg.h>
# include <time.h>
# include <math.h>
# include <float.h>
# include <string.h>
# include <unistd.h>
# include <sys/types.h>
# include <sys/stat.h>

# define PI  M_PI
# define INF 1.0e14

# define EPS 1.0e-14

# define EE 0
# define II 1
# define DD 2
# define dd 3

extern int DEBUG;
extern int max_generations;         // maximum number of generations
extern int number_runs;             // number of experiment runs
extern int run_index;
extern int popsize;                 // population size
extern int number_variable;         // number of variables
extern int number_objective;        // number of objectives
extern int run_index_begin;
extern int run_index_end;
extern double *variable_lowerbound; // variable lower bound
extern double *variable_upperbound; // variable upper bound
extern double eta_c;
extern double eta_m;
extern double pcross_real;
extern double pmut_real;

extern char dummy[50];
extern char problem_name[50];
extern char algorithm_name[50];
extern char analyse_stream[200];

extern int runtime_output;
extern int output_interval;

extern int PF_size;                 // size of the true Pareto-optimal Front
extern double **PF_data;            // true Pareto-optimal front data

extern int analyse_list[100];
enum analyse_name{VAR, FUN, IGD, HV, END};

typedef struct
{
    int rank;
    double *xreal;
    double *obj;
    double crowd_dist;
} individual_real;

typedef struct
{
    individual_real *ind;
} population_real;

typedef struct lists
{
    int index;
    struct lists *parent;
    struct lists *child;
} list;

void insert (list *node, int x);
list* del (list *node);

# endif // Samaritan_GLOBAL_H
