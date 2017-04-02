/*
 * nsga2.c:
 *  This file contains the main procedures of the standard IBEA.
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

# include "../header/metaheuristics.h"
#include "../header/global.h"

#define rho 1.1
#define kappa 0.05
#define indicator 0

double **fitcomp;
int * flag;

int dominates(individual_real* p_ind_a, individual_real* p_ind_b)
{
    int i;
    int a_is_worse = 0;
    int equal = 1;

    for (i = 0; i < number_objective && !a_is_worse; i++)
    {
        a_is_worse = p_ind_a->obj[i] > p_ind_b->obj[i];
        equal = (p_ind_a->obj[i] == p_ind_b->obj[i]) && equal;
    }

    return (!equal && !a_is_worse);
}

double calcHypervolumeIndicator(individual_real* p_ind_a, individual_real* p_ind_b, int d)
/* calculates the hypervolume of that portion of the objective space that
   is dominated by individual a but not by individual b */
{
    double a, b, r, max;
    double volume = 0;

    r = rho * (variable_upperbound[d - 1] - variable_lowerbound[d - 1]);
    max = variable_lowerbound[d - 1] + r;

    print_error(p_ind_a == NULL,1,"error: p_ind_a = NULL\n");
    a = p_ind_a->obj[d - 1];
    if (p_ind_b == NULL)
        b = max;
    else
        b = p_ind_b->obj[d - 1];

    if (d == 1)
    {
        if (a < b)
            volume = (b - a) / r;
        else
            volume = 0;
    }
    else
    {
        if (a < b)
        {
            volume = calcHypervolumeIndicator(p_ind_a, NULL, d - 1) *
                     (b - a) / r;
            volume += calcHypervolumeIndicator(p_ind_a, p_ind_b, d - 1) *
                      (max - b) / r;
        }
        else
        {
            volume = calcHypervolumeIndicator(p_ind_a, p_ind_b, d - 1) *
                     (max - a) / r;
        }
    }

    return (volume);
}


double calcAddEpsIndicator(individual_real* p_ind_a, individual_real* p_ind_b)
/* calculates the maximum epsilon value by which individual a must be
   decreased in all objectives such that individual b is weakly dominated */
{
    int i;
    double r;
    double eps = 0;

    r = variable_upperbound[0] - variable_lowerbound[0];
    eps = (p_ind_a->obj[0] - variable_lowerbound[0]) / r -
          (p_ind_b->obj[0] - variable_lowerbound[0]) / r;
    for (i = 1; i < number_objective; i++)
    {
        double temp_eps;

        r = variable_upperbound[i] - variable_lowerbound[i];
        temp_eps = (p_ind_a->obj[i] - variable_lowerbound[i]) / r -
                   (p_ind_b->obj[i] - variable_lowerbound[i]) / r;
        if (temp_eps > eps)
            eps = temp_eps;
    }

    return (eps);
}


double calcIndicatorValue(individual_real* p_ind_a, individual_real* p_ind_b)
{
    double indicatorValue;

    if (indicator == 0)
        indicatorValue = calcAddEpsIndicator(p_ind_a, p_ind_b);
    else
    {
        if (dominates(p_ind_a, p_ind_b))
            indicatorValue = -calcHypervolumeIndicator(p_ind_a, p_ind_b, number_objective);
        else
            indicatorValue = calcHypervolumeIndicator(p_ind_b, p_ind_a, number_objective);
    }

    return (indicatorValue);
}



void calcFitnessComponents(void * ptr, int size)
{
    population_real * pop = ptr;
    double maxAbsIndicatorValue = 0;
    int i, j;

    /* determine indicator values and their maximum */
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            fitcomp[i][j] = calcIndicatorValue(&pop->ind[i],
                                               &pop->ind[j]);
            if (maxAbsIndicatorValue < fabs(fitcomp[i][j]))
                maxAbsIndicatorValue = fabs(fitcomp[i][j]);
        }
    }

    /* calculate for each pair of invidiuals the corresponding fitness
       component */
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++) {
            fitcomp[i][j] = exp((-fitcomp[i][j]/maxAbsIndicatorValue)/kappa);
        }
    }

    return;
}



void cal_fitnesses(void* ptr, int size)
{
    int i, j;
    double sum;
    population_real * pop = (population_real*) ptr;

    for (i = 0; i < size; i++)
    {
        sum = 0;
        for (j = 0; j < size; j++)
            if (i != j)
                sum += fitcomp[j][i];
        pop->ind[i].crowd_dist = sum;
        //pp_all->ind_array[i]->fitness = sum;
    }

    return;
}



void environmental_selection(void * mixed_ptr,void * new_ptr, int size)
{

    int i, j, worst;
    int new_size = 0;
    population_real * pop = mixed_ptr;
    population_real * new_pop = new_ptr;
    for(i = 0;i< size;i++)
        flag[i] = 0;

    for (i = size - popsize; i > 0; i--)
    {
        for (j = 0; j < size && flag[j] == 1; j++);

        worst = j;

        for (j = j + 1; j < size; j++)
        {
            if (flag[j] != 1)
            {
                if (pop->ind[j].crowd_dist >
                        pop->ind[worst].crowd_dist)
                    worst = j;
            }
        }

        for (j = 0; j < size; j++)
            if (flag[j] != 1)
                pop->ind[j].crowd_dist -= fitcomp[worst][j];

        flag[worst] = 1;
    }

    /* Move remaining individuals to top of array in 'pp_all' */
    for (i = 0; i < size; i++)
    {
        if (flag[i] != 1)
        {
            copy_ind(&pop->ind[i],&new_pop->ind[new_size]);

            new_size++;
        }
    }

    return;
}



void ibea_selection(void * mixed_pop, void *new_pop)
{
    int i;
    int size = 2 * popsize;

    calcFitnessComponents(mixed_pop,size);

    cal_fitnesses(mixed_pop,size);

    environmental_selection(mixed_pop,new_pop,size);

    return;
}


void IBEA (population_real* parent_pop, population_real* offspring_pop, population_real* mixed_pop)
{
    int i;
    int generation;

    generation       = 1;
    evaluation_count = 0;
    printf ("Progress: 1%%");

    // initialize population
    initialize_population_real (parent_pop);

    // population evaluations
    evaluate_population (parent_pop);

    flag = (int *) malloc(2*popsize * sizeof(int));
    fitcomp = (double**) malloc(2*popsize * sizeof(double*));
    for (i = 0; i < 2 * popsize; i++)
        fitcomp[i] = (double*) malloc(2*popsize * sizeof(double*));

    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, generation, 0);
    while (evaluation_count<max_evaluation)
    {
        generation++;
        print_progress ();

        // reproduction (crossover and mutation)
        crossover_real (parent_pop, offspring_pop);

        mutation_real (offspring_pop);

        // population evaluations
        evaluate_population (offspring_pop);

        // environmental selection
        merge (parent_pop, offspring_pop, mixed_pop);

        ibea_selection (mixed_pop, parent_pop);

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, generation, evaluation_count >= max_evaluation);
    }
    for( i = 0;i < 2 * popsize; i++)
        free(fitcomp[i]);
    free(fitcomp);
    free(flag);
    return;
}
