/*
 * spea2.c:
 *  This file contains the main procedures of the standard SPEA2.
 *
 * Authors:
 *  Qi Xu <qixu.student@gmail.com>
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *  Ke Li <k.li@exeter.ac.uk>
 *
 * Institution:
 *  Computational Optimization and Data Analytics (CODA) Group @ University of Exeter
 *
 * Copyright (c) 2017 Qi Xu, Renzhi Chen, Ke Li
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

# include "../header/metaheuristics.h"
#include "../header/global.h"

void SPEA2(population_real *parent_pop, population_real *archive, population_real *mixed_pop) {
    // TODO: implement the framework of SPEA2.
    // TODO: How to deal with the archive?
    // Here, the offspring is actually used as archive, so the name is changed for convenience.
    int archive_size = popsize;
    int total_size = popsize + archive_size;
    int k_neighbor = (int)floor(sqrt((double)total_size));
    int i,j;
    int generation;

    generation = 1;
    evaluation_count = 0;
    printf("Progress: 1%%");

    // initialize population
    initialize_population_real(parent_pop);
    evaluate_population(parent_pop);
    for (i = 0; i < popsize; i++) {
        copy_ind(&(parent_pop->ind[i]), &(archive->ind[i]));
    }

    // track the current evolutionary progress, including population and metrics
    track_evolution(archive, generation, 0);
    while (evaluation_count < max_evaluation) {
        generation++;
        print_progress();

        // reproduction (crossover and mutation)
        crossover_spea2(archive, parent_pop); // Mating selection included here.
        mutation_real(parent_pop);
        evaluate_population(parent_pop);

        // environmental selection
        merge(parent_pop, archive, mixed_pop);
        fitness_spea2(mixed_pop,total_size);
        selection_spea2(mixed_pop,total_size,archive,archive_size);

        // track the current evolutionary progress, including population and metrics
        track_evolution(archive, generation, evaluation_count >= max_evaluation);
    }

    return;
}















