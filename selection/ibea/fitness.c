/*
 * fitness.c:
 *  This file contains the functions for calculating the fitness values for IBEA.
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

# include "../../header/selection.h"

# define rho 1.1
# define kappa 0.05
# define indicator 0

int dominates (individual_real *ind1, individual_real *ind2)
{
    int i;
    int a_is_worse = 0;
    int equal = 1;

    for (i = 0; i < number_objective && !a_is_worse; i++)
    {
        a_is_worse = ind1->obj[i] > ind2->obj[i];
        equal = (ind1->obj[i] == ind2->obj[i]) && equal;
    }

    return (!equal && !a_is_worse);
}

/* calculates the hypervolume of that portion of the objective space that
   is dominated by individual a but not by individual b */
double calcHypervolumeIndicator (individual_real* ind1, individual_real* ind2, int d)
{
    double a, b, r, max;
    double volume = 0;

    r   = rho * (variable_upperbound[d - 1] - variable_lowerbound[d - 1]);
    max = variable_lowerbound[d - 1] + r;

    print_error (ind1 == NULL, 1, "error: ind1 = NULL\n");
    a = ind1->obj[d - 1];
    if (ind2 == NULL)
        b = max;
    else
        b = ind2->obj[d - 1];

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
            volume = calcHypervolumeIndicator(ind1, NULL, d - 1) * (b - a) / r;
            volume += calcHypervolumeIndicator(ind1, ind2, d - 1) * (max - b) / r;
        }
        else
        {
            volume = calcHypervolumeIndicator(ind1, ind2, d - 1) * (max - a) / r;
        }
    }

    return volume;
}

/* calculates the maximum epsilon value by which individual a must be
   decreased in all objectives such that individual b is weakly dominated */
double calcAddEpsIndicator (individual_real *ind1, individual_real *ind2)
{
    int i;
    double r;
    double eps, temp_eps;

    r   = variable_upperbound[0] - variable_lowerbound[0];
    eps = (ind1->obj[0] - variable_lowerbound[0]) / r - (ind2->obj[0] - variable_lowerbound[0]) / r;
    for (i = 1; i < number_objective; i++)
    {
        r = variable_upperbound[i] - variable_lowerbound[i];
        temp_eps = (ind1->obj[i] - variable_lowerbound[i]) / r - (ind2->obj[i] - variable_lowerbound[i]) / r;
        if (temp_eps > eps)
            eps = temp_eps;
    }

    return eps;
}

double calcIndicatorValue (individual_real* ind1, individual_real* ind2)
{
    double indicatorValue;

    if (indicator == 0)
        indicatorValue = calcAddEpsIndicator (ind1, ind2);
    else
    {
        if (dominates (ind1, ind2))
            indicatorValue = -calcHypervolumeIndicator (ind1, ind2, number_objective);
        else
            indicatorValue = calcHypervolumeIndicator (ind2, ind1, number_objective);
    }

    return indicatorValue;
}

void calcFitnessComponents (void *ptr, double **fitcomp, int size)
{
    int i, j;
    double maxAbsIndicatorValue;
    population_real *pop = ptr;

    /* determine indicator values and their maximum */
    maxAbsIndicatorValue = 0;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            fitcomp[i][j] = calcIndicatorValue (&pop->ind[i], &pop->ind[j]);
            if (maxAbsIndicatorValue < fabs (fitcomp[i][j]))
                maxAbsIndicatorValue = fabs (fitcomp[i][j]);
        }
    }

    /* calculate for each pair of individuals the corresponding fitness component */
    for (i = 0; i < size; i++)
        for (j = 0; j < size; j++)
            fitcomp[i][j] = exp ((-fitcomp[i][j] / maxAbsIndicatorValue) / kappa);

    return;
}

/* Assign the fitness values to the solutions within a population */
void cal_fitnesses (void *ptr, double **fitcomp, int size)
{
    int i, j;
    double sum;
    population_real *pop = (population_real*) ptr;

    for (i = 0; i < size; i++)
    {
        sum = 0;
        for (j = 0; j < size; j++)
            if (i != j)
                sum += fitcomp[j][i];
        pop->ind[i].fitness = sum;
    }

    return;
}
