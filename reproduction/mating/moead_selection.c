/*
 * crossover_nsga2.c:
 *  This file contains the functions to perform crossover operations in NSGA-II.
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
#include "../../header/global.h"
#include "../../header/int_vector.h"
#include "../../header/rand.h"
#include "../../header/reproduction.h"
void parent_selection(population_real *parent_pop,individual_real*** parents,int sub_problem_id,int neighbor_type)
{
    struct int_vector *mating_pool;
    int parentsize=3;
    int number_to_select = parentsize-1 ;
    int selected_count;
    int *selected_flag;
    int random;
    int i,j;
    int selected_solution;

    (*parents) = malloc(parentsize * sizeof(individual_real*));
    selected_flag = malloc(popsize * sizeof(int));
    for(i = 0 ; i < popsize ; i++)
        selected_flag[i] = 0;
    selected_count = 0;

    while (selected_count < number_to_select)
    {
        if (neighbor_type == NEIGHBOR)
        {
            random = rnd(0, neighbor_size - 1);
            selected_solution = neighborhood[sub_problem_id][random];
        } else
        {
            selected_solution = rnd(0, popsize - 1);
        }


        if (selected_flag[selected_solution] != 1)
        {
            (*parents)[selected_count] = &(parent_pop->ind[selected_solution]);
            selected_count ++;
            selected_flag[selected_solution] = 1;
        }
    }
    (*parents)[2] =&(parent_pop->ind[sub_problem_id]);

    free(selected_flag);

}
