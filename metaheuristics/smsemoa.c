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
# include "../header/initialization.h"
void crossover_real_sms (population_real *parent_pop, individual_real* offspring1,individual_real* offspring2)
{
    int i;
    int temp;
    int rand;
    int *a1, *a2;

    individual_real *parent1, *parent2;

    a1 = (int *)malloc(popsize * sizeof(int));
    a2 = (int *)malloc(popsize * sizeof(int));
    for (i = 0; i < popsize; i++)
        a1[i] = a2[i] = i;

    for (i = 0; i < popsize; i++)
    {
        rand     = rnd (i, popsize - 1);
        temp     = a1[rand];
        a1[rand] = a1[i];
        a1[i]    = temp;
        rand     = rnd (i, popsize - 1);
        temp     = a2[rand];
        a2[rand] = a2[i];
        a2[i]    = temp;
    }

    parent1 = tournament (&parent_pop->ind[a1[0]], &parent_pop->ind[a1[1]]);
    parent2 = tournament (&parent_pop->ind[a1[2]], &parent_pop->ind[a1[3]]);
    sbx_crossover (parent1, parent2, offspring1, offspring2);
    /*
    parent1 = tournament (&parent_pop->ind[a2[i]], &parent_pop->ind[a2[i + 1]]);
    parent2 = tournament (&parent_pop->ind[a2[i + 2]], &parent_pop->ind[a2[i + 3]]);
    sbx_crossover (parent1, parent2, &offspring_pop->ind[i + 2], &offspring_pop->ind[i + 3]);
    */

    free (a1);
    free (a2);

    return;
}


void SMSEMOA (population_real* parent_pop, population_real* offspring_pop, population_real* mixed_pop)
{
    int i,j;
    int generation;
    FILECONTENTS *f = malloc(sizeof(FILECONTENTS));

    generation       = 1;
    evaluation_count = 0;
    printf ("Progress: 1%%");

    // initialize population
    initialize_population_real (parent_pop);

    // population evaluations
    evaluate_population (parent_pop);

    initialize_nadirpoint(parent_pop);

    // malloc
    i_maxn = number_objective;
    i_maxm =  popsize+1;
    int maxdepth = i_maxn - 2;
    i_fs = malloc(sizeof(FRONT) * maxdepth);
    for ( i = 0; i < maxdepth; i++) {
        i_fs[i].points = malloc(sizeof(POINT) * i_maxm);
        for ( j = 0; j < i_maxm; j++) {
            i_fs[i].points[j].objectives = malloc(sizeof(OBJECTIVE) * (i_maxn - i - 1));
        }
    }
    partial = malloc(sizeof(double) * i_maxm);
    heap = malloc(sizeof(int) * i_maxm);
    stacksize = malloc(sizeof(int) * i_maxm);
    stacks = malloc(sizeof(SLICE*) * i_maxm);
    int maxStackSize = MIN(i_maxn-2,i_slicingDepth(i_maxn))+1;
    for ( i=0; i<i_maxm; i++) {
        stacks[i] = malloc(sizeof(SLICE) * maxStackSize);
        for ( j=1; j<maxStackSize; j++) {
            stacks[i][j].front.points = malloc(sizeof(POINT) * i_maxm);
        }
    }

    fsorted = malloc(sizeof(FRONT) * i_maxn);
    for ( i=0; i<i_maxn; i++) {
        fsorted[i].points = malloc(sizeof(POINT) * i_maxm);
    }
    torder = malloc(sizeof(int*) * MAX(i_maxm,i_maxn));
    tcompare = malloc(sizeof(int*) * i_maxm);
    for ( i=0; i<MAX(i_maxn,i_maxm); i++) {
        torder[i] = malloc(sizeof(int) * i_maxn);
    }
    for ( i=0; i<i_maxm; i++) {
        tcompare[i] = malloc(sizeof(int) * i_maxn);
    }



    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, generation, 0);
    while (evaluation_count<max_evaluation)
    {
        generation ++;
        print_progress ();
        for(j=0;j<popsize;j++) {
            fflush(stdout);
            // reproduction (crossover and mutation)
            crossover_real_sms(parent_pop, &(offspring_pop->ind[0]), &(offspring_pop->ind[1]));

            mutation_ind(&(offspring_pop->ind[0]));

            // population evaluations
            evaluate_individual(&(offspring_pop->ind[0]));


            update_nadir_point(&(offspring_pop->ind[0]));

            int kk;
            //for(kk=0;kk<number_objective;kk++)
            //    printf("%lf ",nadir_point[kk]);
            //printf("\n");

            // environmental selection
            merge(parent_pop, offspring_pop, mixed_pop);

            //fill_nondominated_sort (parent_pop, mixed_pop);
            fill_hv_sort(f, parent_pop, mixed_pop, popsize + 1);
        }
            // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, generation, evaluation_count >= max_evaluation);
    }

    for ( i = 0; i < maxdepth; i++) {
        for ( j = 0; j < i_maxm; j++) {
            free(i_fs[i].points[j].objectives);
        }
        free(i_fs[i].points);

    }
    free(i_fs);

    for ( i=0; i<i_maxm; i++) {
        free(stacks[i]);
    }
    free(partial);
    free(heap);
    free(stacksize);
    free(stacks);

    for ( i=0; i<i_maxn; i++) {
        free(fsorted[i].points);
    }
    free(fsorted);
    for ( i=0; i<MAX(i_maxn,i_maxm); i++) {
        free(torder[i]);
    }
    for ( i=0; i<i_maxm; i++) {
        free(tcompare[i]);
    }

    free(torder);
    free(tcompare);

    return;
}
