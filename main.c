/*
 * main.c:
 *  This is the main procedures of a general EMO algorithm (generational evolution model).
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

# include "header/global.h"
# include "header/rand.h"
# include "header/metaheuristics.h"
# include "header/population.h"
# include "header/memory.h"
# include "header/analyse.h"

int DEBUG;

int number_runs;             // number of experiment runs
int popsize;                 // population size
int number_variable;         // number of variables
int number_objective;        // number of objectives
int run_index;               // index of the current experiment
int run_index_begin;
int run_index_end;
double *variable_lowerbound; // variable lower bound
double *variable_upperbound; // variable upper bound
double eta_c;
double eta_m;
double pcross_real;
double pmut_real;
int neighbor_size;
int reference_size;
double neighborhood_selection_probability;
double DEFAULT_CR;
double DEFAULT_F;
double DEFAULT_K;

double ** lambda;
int **neighborhood;
double *ideal_point;
int * permutation;

char dummy[BUFSIZE_S];
char problem_name[BUFSIZE_S];
char algorithm_name[BUFSIZE_S];
char analyse_stream[BUFSIZE_L];

int runtime_output;
int output_interval;
int evaluation_count;
int max_evaluation;
int PF_size;                 // size of the true Pareto-optimal Front
double **PF_data;            // true Pareto-optimal front data

int analyse_list[BUFSIZE_S];

int main(int argc, char *argv[])
{
    int i;

    // initialize parameter settings
    init_real ("config.txt");

    population_real* parent_pop;
    population_real* offspring_pop;
    population_real* mixed_pop;
    parent_pop    = (population_real *) malloc (sizeof(population_real));
    offspring_pop = (population_real *) malloc (sizeof(population_real));
    mixed_pop     = (population_real *) malloc (sizeof(population_real));
    allocate_memory_pop (parent_pop, popsize);
    allocate_memory_pop (offspring_pop, popsize);
    allocate_memory_pop (mixed_pop, 2 * popsize);

    randomize ();

    // run experiments
    for (run_index = run_index_begin; run_index <= run_index_end; run_index++) {
        printf ("The %d run ...\n", run_index);
        if (!strcmp(algorithm_name, "NSGA2"))NSGA2 (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp(algorithm_name, "MOEAD"))MOEAD (parent_pop, offspring_pop, mixed_pop);

        else
        {
            print_error (1, 2, "UNKNOWN algorithm:", algorithm_name);
        }
    }
    //analyse_all ();

    // free memory
    if (number_variable != 0)
    {
        free (variable_lowerbound);
        free (variable_upperbound);
    }
    deallocate_memory_pop (parent_pop, popsize);
    deallocate_memory_pop (offspring_pop, popsize);
    deallocate_memory_pop (mixed_pop, 2 * popsize);
    free (parent_pop);
    free (offspring_pop);
    free (mixed_pop);

    for (i = 0; i < PF_size; i++)
        free (PF_data[i]);
    free (PF_data);

    return 0;
}
