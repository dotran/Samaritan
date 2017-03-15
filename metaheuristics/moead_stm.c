/*
 * moead_stm.c:
 *  This file contains the all procedures of the MOEAD-STM.
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
# include "../header/problems.h"
# include "../header/analyse.h"
# include "../header/initialization.h"

extern int* idx ;
extern double* nicheCount;
extern struct double_with_index** solMatrix;
extern double** distMatrix;
extern double** fitnessMatrix;
extern struct double_with_index** subpMatrix;
extern int* statusWoman;
extern int* next;


void stm_selection(population_real* parent_pop,population_real* mixed_pop)
{

    int i,j;

    // Calculate the preference values of solution matrix
    for (i = 0; i < 2*number_weight; i++)
    {
        int minIndex = 0;
        for (j = 0; j < number_weight; j++) {
            fitnessMatrix[i][j] = fitnessFunction(&(mixed_pop->ind[i]), lambda[j]);
            distMatrix[i][j]  	= calculateDistance2(&(mixed_pop->ind[i]), lambda[j]);
            if (distMatrix[i][j] < distMatrix[i][minIndex])
                minIndex = j;
        }
        nicheCount[minIndex] = nicheCount[minIndex] + 1;
    }

    // calculate the preference values of subproblem matrix and solution matrix
    for (i = 0; i < 2*number_weight; i++)
    {
        for (j = 0; j < number_weight; j++) {
            subpMatrix[j][i].x = fitnessFunction(&(mixed_pop->ind[i]), lambda[j]);
            subpMatrix[j][i].idx = i;
            solMatrix[i][j].x = distMatrix[i][j] + nicheCount[j];
            solMatrix[i][j].idx = j;
        }
    }

    for ( i = 0 ; i < number_weight ; i++)
    {
        qsort(subpMatrix[i],2*number_weight,sizeof(struct double_with_index),double_with_index_greater_cmp);
    }
    for ( i = 0 ; i < 2 * number_weight ; i++)
    {
        qsort(solMatrix[i],number_weight,sizeof(struct double_with_index),double_with_index_greater_cmp);
    }
    /*
    printf("subpMatrix:\n");
    for(i =0;i<popsize;i++) {
        for (j = 0; j < 2 * popsize; j++)
            printf("%lf(%d) ", subpMatrix[i][j].x, subpMatrix[i][j].idx);
        printf("\n");
    }
    printf("solMatrix:\n");
    for(i =0;i<2*popsize;i++) {
        for (j = 0; j < popsize; j++)
            printf("%lf(%d) ", solMatrix[i][j].x, solMatrix[i][j].idx);
        printf("\n");
    }
    */
    stableMatching(idx,statusWoman,next,subpMatrix, solMatrix, number_weight, number_weight*2);

    for (i = 0; i < number_weight; i++)
    {

        copy_ind(&(mixed_pop->ind[idx[i]]), &(parent_pop->ind[i]));
    }
}


void MOEAD_STM (population_real* parent_pop, population_real* offspring_pop, population_real* mixed_pop)
{

    int i;
    int generation;
    int subproblem_id, neighbor_type;

    initialize_uniform_weight ();

    idx = malloc(sizeof(int)*number_weight);
    nicheCount = malloc(sizeof(double)*number_weight);
    distMatrix    = malloc(sizeof(double*)*number_weight*2);
    fitnessMatrix = malloc(sizeof(double*)*number_weight*2);
    statusWoman = malloc(sizeof(int)*number_weight*2);
    next = malloc(sizeof(int)*number_weight*2);
    solMatrix = malloc(sizeof(struct double_with_index*)*number_weight*2);
    subpMatrix = malloc(sizeof(struct double_with_index*)*number_weight);

    for (i = 0; i < 2*number_weight; i++)
    {
        solMatrix[i] = malloc(sizeof(struct double_with_index)*number_weight);
        distMatrix[i]    = malloc(sizeof(double)*number_weight);
        fitnessMatrix[i] = malloc(sizeof(double)*number_weight);
    }

    for ( i = 0; i < number_weight; i++)
    {
        subpMatrix[i] = malloc(sizeof(struct double_with_index)*number_weight*2);;
    }

    generation       = 1;
    evaluation_count = 0;
    printf ("Progress: 1%%");


    // initialize population
    initialize_population_real (parent_pop);

    // population evaluations
    evaluate_population (parent_pop);

    initialize_neighborhood ();

    initialize_idealpoint (parent_pop);

    initialize_nadirpoint(parent_pop);

    permutation = malloc (number_weight * sizeof(int));

    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, generation, 0);

    while (evaluation_count < max_evaluation)
    {
        print_progress (generation);

        random_permutation (permutation,number_weight);
        for (i = 0; i < number_weight; i++)
        {
            subproblem_id = permutation[i];
            crossover_moead_real (parent_pop, &(offspring_pop->ind[i]), subproblem_id, &neighbor_type);
            mutation_ind (&(offspring_pop->ind[i]));
            evaluate_individual (&(offspring_pop->ind[i]));
            update_ideal_point (&(offspring_pop->ind[i]));
            update_nadir_point(&(offspring_pop->ind[i]));
        }

        merge (parent_pop, offspring_pop, mixed_pop);

        stm_selection(parent_pop,mixed_pop);

        generation++;

        track_evolution (parent_pop, generation, evaluation_count >= max_evaluation);

    }

    // free all malloc memory
    for (i = 0; i < 2*number_weight; i++)
    {
        free(solMatrix[i]);
        free(distMatrix[i]);
        free(fitnessMatrix[i]);
    }

    for (i = 0; i < number_weight; i++)
    {
        free(subpMatrix[i]);
    }
    free(permutation);
    free(idx) ;
    free(nicheCount);
    free(solMatrix);
    free(distMatrix);
    free(fitnessMatrix);
    free(statusWoman);
    free(next);
    free (permutation);
    moead_free ();
    return;
}