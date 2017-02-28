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

void NSGA2 (population_real* parent_pop, population_real* offspring_pop, population_real* mixed_pop)
{
    int i;
    int generation;

    generation = 1;
    printf ("Progress: 1%%");

    // initialize population
    initialize_population_real (parent_pop);

    evaluation_count = 0;

    // population evaluations
    evaluate_population (parent_pop);

    double t[10];
    for(i=0;i<10;i++)t[i] =0;
    clock_t current;

    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, generation, 0);
    while(evaluation_count<max_evaluation)
    {

        generation ++;
        print_progress (generation);

        // reproduction (crossover and mutation)
        current = clock();
        crossover_real (parent_pop, offspring_pop);
        t[0] += 1.0 * (clock()-current)/CLOCKS_PER_SEC;
        current = clock();
        mutation_real (offspring_pop);
        t[1] += 1.0 * (clock()-current)/CLOCKS_PER_SEC;
        current = clock();
        // population evaluations
        evaluate_population (offspring_pop);
        t[2] += 1.0 * (clock()-current)/CLOCKS_PER_SEC;
        current = clock();
        // environmental selection
        merge (parent_pop, offspring_pop, mixed_pop);
        fill_nondominated_sort (parent_pop, mixed_pop);
        t[3] += 1.0 * (clock()-current)/CLOCKS_PER_SEC;

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, generation,evaluation_count>=max_evaluation);



    }
    printf("%lf,%lf,%lf,%lf\n",t[0],t[1],t[2],t[3]);

    return;
}
