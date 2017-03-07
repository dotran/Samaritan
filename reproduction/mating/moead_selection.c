/*
 * moead_selection.c:
 *  This file contains the selecting functions for moead.
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
#include "../../header/rand.h"
#include "../../header/reproduction.h"
#include "../../header/vector.h"

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



int* tour_selection(int depth) {

    // selection based on utility
    int i,k,n;
    int selected_size;
    int candidate_size;

    for (k = 0; k < number_objective; k++)
        int_vector_pushback(selected,k);
    selected_size = int_vector_size(selected);

    for (n = number_objective; n < popsize; n++)
        int_vector_pushback(candidate,n);
    candidate_size = int_vector_size(candidate);

    while (selected_size < (int) (popsize / 5.0)) {

        int best_idd = (int) (rndreal(0,1) * candidate_size);
        int i2;
        int best_sub = int_vector_get(candidate,best_idd+1);
        int s2;
        for (i = 1; i < depth; i++) {
            i2 = (int) (rndreal(0,1) * candidate_size);
            s2 = int_vector_get(candidate,i2+1);
            // System.out.println("Candidate: "+i2);
            if (utility[s2] > utility[best_sub]) {
                best_idd = i2;
                best_sub = s2;
            }
        }
        int_vector_pushback(selected, best_sub);
        selected_size ++;

        int_vector_remove(candidate,best_idd+1);
        candidate_size --;

    }
}

void comp_utility(population_real* pop,population_real* saved_values) {
    double f1, f2, uti, delta;
    int n;
    for (n = 0; n < popsize; n++) {
        f1 = fitnessFunction(&(pop->ind[n]), lambda[n]);
        f2 = fitnessFunction(&(saved_values->ind[n]), lambda[n]);
        delta = f2 - f1;
        if (delta > 0.001)
            utility[n] = 1.0;
        else {
            // uti = 0.95*(1.0+delta/0.001)*utility_[n];
            uti = (0.95 + (0.05 * delta / 0.001)) * utility[n];
            utility[n] = uti < 1.0 ? uti : 1.0;
        }
        copy_ind(&(pop->ind[n]), &(saved_values->ind[n]));
    }

}