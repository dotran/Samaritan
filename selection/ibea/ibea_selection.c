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
#include "../../header/global.h"

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
    int size;

    size = 2 * popsize;
    calcFitnessComponents (mixed_pop, fitcomp, size);
    cal_fitnesses (mixed_pop, fitcomp, size);
    environmental_selection (mixed_pop, new_pop, flag, fitcomp, size);

}

void pbea_selection(void *mixed_pop, void *new_pop, int *flag, double **fitcomp, double* weights, double* reference_point, double specificity) {
    int size = 2 * popsize;
    pbea_calcFitnessComponents (mixed_pop, fitcomp, size, weights, reference_point, specificity);
    cal_fitnesses (mixed_pop, fitcomp, size);
    environmental_selection (mixed_pop, new_pop, flag, fitcomp, size);
}

void cibea_selection (void *mixed_ptr, void *new_ptr, int *flag, double **fitcomp) {
    int i;
    int size;
    int cv_size = 0;
    size = 2 * popsize;
    population_real * mixed_pop = mixed_ptr;
    population_real * new_pop = new_ptr;
    for (i = 0; i < size; i++)
    {
        if(mixed_pop->ind[i].cv>-EPS)
            cv_size++;

    }

    if(cv_size>popsize)
    {
        int left = 0;
        int right = size -1;
        individual_real temp;

        while(left<right)
        {
            while(mixed_pop->ind[left].cv>-EPS&&left<right)
                left++;
            while(mixed_pop->ind[right].cv<-EPS&&left<right)
                right--;
            if(left<right)
            {
                //printf("[%d]%lf <- [%d]%lf\n",left,mixed_pop->ind[left].cv,right,mixed_pop->ind[right].cv);
                fflush(stdout);
                copy_ind(&mixed_pop->ind[right], &mixed_pop->ind[left]);
                right --;
                left ++;
            }
        }

        calcFitnessComponents (mixed_ptr, fitcomp, cv_size);
        cal_fitnesses (mixed_ptr, fitcomp, cv_size);
        environmental_selection (mixed_ptr, new_pop, flag, fitcomp, cv_size);
    }
    else
    {
        int selected = 0;
        struct double_with_index * temp;
        temp = malloc(sizeof(struct double_with_index) * size);
        for ( i = 0; i < size; i++)
        {
            if(mixed_pop->ind[i].cv>-EPS)
            {
                copy_ind(&mixed_pop->ind[i], &new_pop->ind[selected]);
                selected++;
            }
            else
            {
                temp[i-selected].idx = i;
                temp[i-selected].x = mixed_pop->ind[i].cv;
            }
        }
        qsort (temp, size- cv_size, sizeof(struct double_with_index), double_with_index_smaller_cmp);
        for(i=0;i<popsize-cv_size;i++)
        {
            copy_ind(&mixed_pop->ind[temp[i].idx],&new_pop->ind[selected]);
            selected++;
        }
        free(temp);

    }

    return;
}
