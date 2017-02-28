/*
 * nsga2.c:
 *  This file contains the main procedures of the standard NSGA-II.
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
# include "../header/utility.h"
# include "./util/moeadutil.h"

void MOEAD (population_real* pop, population_real* offspring_pop, population_real* mixed_pop)
{
    int i,j;
    int generation;
    individual_real * offspring =&( offspring_pop->ind[0]);
    generation = 1;
    printf ("Progress: 1%%");

    // initialize population

    initialize_population_real (pop);

    initialize_uniform_weight ();

    initialize_neighborhood ();

    initialize_idealpoint (pop);

    evaluation_count = 0;

    evaluate_population (pop);

    track_evolution (pop, generation, 0);

    while(evaluation_count<max_evaluation) {

        generation ++;
        printf("\n %d \n",generation);
        print_progress (generation);

        permutation = malloc(popsize*sizeof(int));
        random_permutation(permutation,popsize);

        for(i = 0 ; i < popsize; i++)
        {
            int neighbor_type;
            int sub_problem_id = permutation[i];

            crossover_moead_real (pop, offspring,sub_problem_id,&neighbor_type);

            mutation_ind (offspring);

            evaluate_individual (offspring);

            update_ideal_point(offspring);

            update_neighborhood(pop,offspring, sub_problem_id, neighbor_type);

            track_evolution (pop, generation,evaluation_count>=max_evaluation);
        }
        free(permutation);
    }
    moead_free();
/*
    evaluation_count = 0;

    // population evaluations
    evaluate_population (parent_pop);

    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, generation, 0);
    while(evaluation_count<max_evaluation)
    {

        generation ++;
        print_progress (generation);

        // reproduction (crossover and mutation)
        crossover_real (parent_pop, offspring_pop);
        mutation_real (offspring_pop);

        // population evaluations
        evaluate_population (offspring_pop);

        // environmental selection
        merge (parent_pop, offspring_pop, mixed_pop);
        fill_nondominated_sort (parent_pop, mixed_pop);


        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, generation,evaluation_count>=max_evaluation);



    }
*/
    return;
}

