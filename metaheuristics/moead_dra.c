/*
 * moead_dra.c:
 *  This file contains the main procedures of the standard MOEAD.
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
# include "../header/global.h"
# include "../header/population.h"
# include "../header/reproduction.h"
# include "../header/selection.h"
# include "../header/problems.h"
# include "../header/analyse.h"
# include "../header/memory.h"

void MOEAD_DRA (population_real* pop, population_real* offspring_pop, population_real* mixed_pop)
{
    int i, j;
    int generation;
    int subproblem_id, neighbor_type;
    population_real* saved_values;

    generation       = 1;
    evaluation_count = 0;
    printf ("Progress: 1%%");

    // initialization process
    initialize_uniform_weight ();
    initialize_neighborhood ();
    initialize_population_real (pop);
    evaluate_population (pop);
    initialize_idealpoint (pop);

    track_evolution (pop, generation, 0);

    saved_values = malloc (popsize * sizeof(population_real));
    allocate_memory_pop (saved_values, popsize);

    utility = malloc (popsize * sizeof(double));
    frequency = malloc (popsize * sizeof(int));

    for (i = 0; i < popsize; i++) {
        utility[i] = 1.0;
        frequency[i] = 0;
    }

    permutation = malloc (popsize * sizeof(int));
    individual_real* offspring = &(offspring_pop->ind[0]);

    while (evaluation_count < max_evaluation) {

        //create empty head for selected and candidate
        selected = malloc(sizeof(struct int_vector));
        candidate = malloc(sizeof(struct int_vector));
        selected->value = INT_MIN;
        selected->next=NULL;
        candidate->value = INT_MIN;
        candidate->next=NULL;

        print_progress (generation);
        random_permutation (permutation,popsize);
        tour_selection(neighbor_size);
        for (i = 0; i < int_vector_size(selected); i++)
        {

            int n = int_vector_get(selected,i+1);
            frequency[n]++;

            subproblem_id = permutation[i];

            crossover_moead_real (pop, offspring, n, &neighbor_type);
            mutation_ind (offspring);
            evaluate_individual (offspring);

            update_ideal_point (offspring);

            update_neighborhood (pop,offspring, n, neighbor_type);
        }

        generation ++;
        if (generation % 30 == 0) {
            comp_utility(pop,saved_values);
        }
        track_evolution (pop, generation, evaluation_count >= max_evaluation);

        int_vector_free(selected);
        int_vector_free(candidate);
    }

    free (permutation);
    moead_free ();
    free(utility);
    free(frequency);
    deallocate_memory_pop(saved_values,popsize);
    free(saved_values);

    return;
}
