/*
 * ibea_selection.c:
 *  This file contains the environmental selection function for IBEA.
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

void environmental_selection (void *mixed_ptr, void *new_ptr, int *flag, double **fitcomp, int size)
{
    int i, j;
    int worst, new_size;

    population_real *pop     = mixed_ptr;
    population_real *new_pop = new_ptr;

    for (i = 0; i < size; i++)
        flag[i] = 0;

    for (i = size - popsize; i > 0; i--)
    {
        for (j = 0; j < size && flag[j] == 1; j++);

        worst = j;

        for (j = j + 1; j < size; j++)
        {
            if (flag[j] != 1)
            {
                if (pop->ind[j].fitness >
                    pop->ind[worst].fitness)
                    worst = j;
            }
        }

        for (j = 0; j < size; j++)
            if (flag[j] != 1)
                pop->ind[j].fitness -= fitcomp[worst][j];

        flag[worst] = 1;
    }

    /* Move remaining individuals to top of array in 'pp_all' */
    new_size = 0;
    for (i = 0; i < size; i++)
    {
        if (flag[i] != 1)
        {
            copy_ind (&pop->ind[i], &new_pop->ind[new_size]);
            new_size++;
        }
    }

    return;
}

void ibea_selection (void *mixed_pop, void *new_pop, int *flag, double **fitcomp)
{
    int i;
    int size;

    size = 2 * popsize;
    calcFitnessComponents (mixed_pop, fitcomp, size);
    cal_fitnesses (mixed_pop, fitcomp, size);
    environmental_selection (mixed_pop, new_pop, flag, fitcomp, size);

    return;
}
