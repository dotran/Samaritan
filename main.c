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

# include "time.h"
# include "header/global.h"
# include "header/rand.h"
# include "header/metaheuristics.h"
# include "header/population.h"
# include "header/memory.h"

// global variables
int DEBUG;

int algorithm_index;         // index of the running algorithm
int problem_index;           // index of the test problem
int max_generations;         // maximum number of generations
int number_runs;             // number of experiment runs
int popsize;                 // population size
int number_variable;         // number of variables
int number_objective;        // number of objectives
double *variable_lowerbound; // variable lower bound
double *variable_upperbound; // variable upper bound
double eta_c;
double eta_m;
double pcross_real;
double pmut_real;

int main(int argc, char *argv[])
{
    // initialization
    init_real();

    population_real* parent_pop;
    population_real* offspring_pop;
    population_real* mixed_pop;
    parent_pop    = (population_real *) malloc (sizeof(population_real));
    offspring_pop = (population_real *) malloc (sizeof(population_real));
    mixed_pop     = (population_real *) malloc (sizeof(population_real));
    allocate_memory_pop (parent_pop, popsize);
    allocate_memory_pop (offspring_pop, popsize);
    allocate_memory_pop (mixed_pop, 2 * popsize);

    randomize();

    // run experiments
    switch (algorithm_index)
    {
        case 1:
            NSGA2(parent_pop, offspring_pop, mixed_pop);
            break;
        default:
            printf("Please specify an metaheuristic to run.\n");
            exit(1);
    }

    // performance assessment
    // TODO

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

    return 0;
}
