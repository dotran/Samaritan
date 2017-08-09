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

# include "../../header/selection.h"
#include "../../header/global.h"

double search_k_minimum (double *a, int k, int size)
{
    int i, j;
    int temp_int;
    double k_min, temp;

    print_error(k > size, 1, "Error in search_k_minimum: k must be smaller or equal to size!");

    for (i = 0; i < k; i++)
    {
        k_min    = a[i];
        temp_int = i;

        for (j = i + 1; j < size; j++)
        {
            if (k_min > a[j])
            {
                k_min    = a[j];
                temp_int = j;
            }
        }
        if (i != temp_int)
        {
            temp        = a[i];
            a[i]        = a[temp_int];
            a[temp_int] = temp;
        }
    }

    return (k_min);
}

void fitness_spea2 (population_real *pop, int total_size, int k_min, int *dominated_Num, int *bedominated_Num, int **dominated_Matrix, int *R_i, double **distance_Matrix, double* D_i, double *kth_distance)
{
    int i, j, flag;

    for (i = 0; i < total_size; i++)
    {
        dominated_Num[i]   = 0;
        bedominated_Num[i] = 0;
        R_i[i]             = 0;
        D_i[i]             = 0.0;
    }

    // Calculate R(i)
    for (i = 0; i < total_size - 1; i++)
    {
        for (j = i + 1; j < total_size; j++)
        {
            flag = check_dominance(&(pop->ind[i]), &(pop->ind[j]));
            dominated_Matrix[i][j] = flag;
            dominated_Matrix[j][i] = -(flag);

            if (flag == 1)
            {
                dominated_Num[i]   += 1;
                bedominated_Num[j] += 1;
            }
            if (flag == -1)
            {
                bedominated_Num[i] += 1;
                dominated_Num[j]   += 1;
            }
        }
    }
    for (i = 0; i < total_size; i++)
    {
        for (j = 0; j < total_size; j++)
        {
            if (dominated_Matrix[i][j] == -1)
            {
                R_i[i] += dominated_Num[j];
            }
        }
    }

    // Calculate D(i)
    for (i = 0; i < total_size; i++)
    {
        distance_Matrix[i][i] = INF;
        for (j = i + 1; j < total_size; j++)
        {
            distance_Matrix[i][j] = distance_Matrix[j][i] = euclidian_distance(pop->ind[i].xreal, pop->ind[j].xreal, number_variable);
        }
    }
    for (i = 0; i < total_size; i++)
    {
        kth_distance[i] = search_k_minimum(distance_Matrix[i], k_min, total_size);
    }
    for (i = 0; i < total_size; i++)
    {
        D_i[i] = 1.0 / (2.0 + kth_distance[i]);
    }
    for (i = 0; i < total_size; i++)
    {
        pop->ind[i].fitness = R_i[i] + D_i[i];
    }
}

void selection_spea2 (population_real *mixed_pop,int total_size,population_real *archive,int archive_size, individual_real *temp_ind, population_real *temp_pop, double **distance_Matrix)
{
    int i, j, flag;
    int num_nondominated;
    individual_real *ind, *ind1, *ind2;
    double min, temp;

    num_nondominated = 0;
    for (i = 0; i < total_size; i++) // Count num_nondominated
    {
        ind = &(mixed_pop->ind[i]);
        if (ind->fitness < 1.0)
        {
            num_nondominated += 1;
        }
    }
    if (num_nondominated <= archive_size)
    {
        // Sorting the list should be far more efficient
        // TODO: Implement some kind of benchmarks to measure average running time of implementations

        for (i = 0; i < archive_size; i++)
        {
            ind1 = &(mixed_pop->ind[i]);
            min = ind1->fitness;
            flag = i;

            for (j = i + 1; j < total_size; j++)
            {
                ind2 = &(mixed_pop->ind[j]);
                temp = ind2->fitness;
                if (min > temp)
                {
                    flag = j;
                    min  = temp;
                }
            }
            copy_ind (&(mixed_pop->ind[i]), temp_ind);
            copy_ind (&(mixed_pop->ind[flag]), &(mixed_pop->ind[i]));
            copy_ind (temp_ind, &(mixed_pop->ind[flag]));
            copy_ind (&(mixed_pop->ind[i]), &(archive->ind[i]));
        }
    }
    else // Requires truncation: num_nondominated > archive_size
    {
        list* start = (list *) malloc(sizeof(list));
        start->index = NULL;
        list* end = start;

        int oa_size = 0; // oa - oversized archive

        // Copy non-dominated solutions into a list
        list* non_dom_individual;
        for (i = 0; i < total_size; i++)
        {
            if (mixed_pop->ind[i].fitness < 1.0) {

                copy_ind (&(mixed_pop->ind[i]), &(temp_pop->ind[oa_size]));

                non_dom_individual = (list *) malloc(sizeof(list));
                non_dom_individual->index = i;
                non_dom_individual->index2 = oa_size;
                non_dom_individual->parent = end;
                non_dom_individual->child = NULL;

                end->child = non_dom_individual;

                end = non_dom_individual;

                oa_size ++;
            }
        }

        // Calculate the distances between the non-dominated solutions
        for (i = 0; i < oa_size; i++)
        {
            distance_Matrix[i][i] = INF;
            for (j = i + 1; j < oa_size; j++)
            {
                distance_Matrix[i][j] = distance_Matrix[j][i] = euclidian_distance(temp_pop->ind[i].xreal, temp_pop->ind[j].xreal, number_variable);
            }
        }

        int temp_i = 0;
        double temp_val;

        while(oa_size > archive_size) {
            temp_val = INF;

            // Find the closest element
            for(i = 0; i < oa_size; i++) {
                for(j = 0; j < oa_size; j++) {
                    if(distance_Matrix[i][j] < temp_val) {
                        temp_i = i;
                        temp_val = distance_Matrix[i][j];
                    }
                }
            }

            // Remove it from the list
            list* temp_list = start;
            while(temp_list->child != NULL) {
                temp_list = temp_list->child;

                if(temp_list->index2 == temp_i) {
                    del(temp_list);
                    break;
                }
            }

            oa_size --;
        }

        // Copy from list to archive
        list* temp_list = start;

        for(i = 0; i < oa_size; i++) {
            temp_list = temp_list->child;
            copy_ind (&(mixed_pop->ind[temp_list->index]), &(archive->ind[i]));

            free(temp_list->parent);
        }

        free(temp_list);
    }

}









