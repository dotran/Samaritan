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
#include "../../header/memory.h"
#include "../../header/global.h"

double search_k_minimum(double *a, int k, int size) {
    int i, j;
    int temp_int;
    double k_min, temp;

    if (k > size) {
        printf("函数search_k_minimum出错! 退出!");
    }

    for (i = 0; i < k; i++) {
        k_min = a[i];
        temp_int = i;

        for (j = i + 1; j < size; j++) {
            if (k_min > a[j]) {
                k_min = a[j];
                temp_int = j;
            }
        }
        if (i != temp_int) {
            temp = a[i];
            a[i] = a[temp_int];
            a[temp_int] = temp;
        }
    }
    return (k_min);
}

void fitness_spea2(population_real *pop, int size) {
    int i, j, k_min;
    int flag;

    int *dominated_Num;            // 用于记录种群pop中每一个个体所支配种pop中的个体数
    int *bedominated_Num;        // 用于记录种群pop中每一个个体被pop中个体所支配的数目
    int **dominated_Matrix;        // 矩阵中的每一个元素a[i][j]的取值, 1表示i支配j, -1表示j支配i, 0表示j与i互不支配
    int *R_i;                    // 用于记录种群pop中每一个个体的R(i)

    double **distance_Matrix;    // 矩阵中的每一个元素a[i][j]的取值表示个体i与j
    double d;
    double *D_i;                // 用于记录种群pop中每一个个体的D(i)
    double *kth_distance;
    individual_real *ind;

    dominated_Num = (int *) malloc(size * sizeof(int));
    bedominated_Num = (int *) malloc(size * sizeof(int));
    R_i = (int *) malloc(size * sizeof(int));
    D_i = (double *) malloc(size * sizeof(double));
    kth_distance = (double *) malloc(size * sizeof(double));
    dominated_Matrix = (int **) malloc(size * sizeof(int *));
    distance_Matrix = (double **) malloc(size * sizeof(double *));

    for (i = 0; i < size; i++) {
        dominated_Matrix[i] = (int *) malloc(size * sizeof(int));
        distance_Matrix[i] = (double *) malloc(size * sizeof(double));
    }

    k_min = (int) floor(sqrt(size));    // 第k个最近的个体

    for (i = 0; i < size; i++) {
        dominated_Num[i] = 0;
        bedominated_Num[i] = 0;
        R_i[i] = 0;
        D_i[i] = 0.0;
    }

    // Calculate R(i)
    for (i = 0; i < size - 1; i++) {
        for (j = i + 1; j < size; j++) {
            flag = check_dominance(&(pop->ind[i]), &(pop->ind[j]));
            dominated_Matrix[i][j] = flag;
            dominated_Matrix[j][i] = -(flag);

            if (flag == 1) {
                dominated_Num[i] += 1;
                bedominated_Num[j] += 1;
            }
            if (flag == -1) {
                bedominated_Num[i] += 1;
                dominated_Num[j] += 1;
            }
        }
    }
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            if (dominated_Matrix[i][j] == -1) {
                R_i[i] += dominated_Num[j];
            }
        }
    }

    // Calculate D(i)
    for (i = 0; i < size; i++)
    {
        distance_Matrix[i][i] = INF;
        for (j = i + 1; j < size; j++) {
            d = euclidian_distance(pop->ind[i].xreal, pop->ind[j].xreal, number_variable);
            distance_Matrix[i][j] = d;
            distance_Matrix[j][i] = d;
        }
    }
    for (i = 0; i < size; i++) {
        kth_distance[i] = search_k_minimum(distance_Matrix[i], k_min, size);
    }
    for (i = 0; i < size; i++) {
        D_i[i] = 1.0 / (2.0 + kth_distance[i]);
    }
    for (i = 0; i < size; i++) {
        ind = &(pop->ind[i]);
        ind->fitness = R_i[i] + D_i[i];
    }

    free(dominated_Num);
    free(bedominated_Num);
    free(R_i);
    free(D_i);
    free(kth_distance);
    for (i = 0; i < size; i++) {
        free(dominated_Matrix[i]);
        free(distance_Matrix[i]);
    }
    free(distance_Matrix);
    free(dominated_Matrix);
}

void truncate_pop(population_real *tmppop, int tmppop_size, population_real *pop2, int archive_size) {
    int i, j, flag, flag2;
    list *elist, *temp_list, *plist;
    double **distance_Matrix, d, temp1, temp2;

    elist = (list *) malloc(sizeof(list));
    distance_Matrix = (double **) malloc(tmppop_size * sizeof(double));

    for (i = 0; i < tmppop_size; i++) {
        distance_Matrix[i] = (double *) malloc(tmppop_size * sizeof(double));
    }

    elist->index = -1;
    elist->parent = NULL;
    elist->child = NULL;
    temp_list = elist;

    for (i = 0; i < tmppop_size; i++) {
        insert(temp_list, i);
        temp_list = temp_list->child;
    }

    for (i = 0; i < tmppop_size; i++)
    {
        distance_Matrix[i][i] = INF;

        for (j = i + 1; j < tmppop_size; j++) {
            d = euclidian_distance(tmppop->ind[i].xreal, tmppop->ind[j].xreal, number_variable);
            distance_Matrix[i][j] = d;
            distance_Matrix[j][i] = d;
        }
    }

    for (i = 0; i < tmppop_size - archive_size; i++) {
        temp_list = elist->child;
        temp1 = 0.0;

        while (temp_list != NULL) {
            flag = 0;
            temp1 = distance_Matrix[temp_list->index][0];

            for (j = 1; j < tmppop_size; j++) {
                if (distance_Matrix[temp_list->index][j] < temp1) {
                    temp1 = distance_Matrix[temp_list->index][j];
                    flag = j;
                }
            }
            temp_list->index2 = flag;
            temp_list = temp_list->child;
        }

        temp_list = elist->child;
        temp1 = distance_Matrix[temp_list->index][temp_list->index2];
        temp_list = elist->child;
        plist = temp_list;
        while (temp_list != NULL)
        {
            if (distance_Matrix[temp_list->index][temp_list->index2] < temp1) {
                temp1 = distance_Matrix[temp_list->index][temp_list->index2];
                plist = temp_list;
                temp_list = temp_list->child;
            } else {
                temp_list = temp_list->child;
            }
        }

        if (plist != NULL) {
            flag2 = plist->index;
            if (plist->index == 0) {
                temp1 = distance_Matrix[plist->index2][plist->index + 1];

                for (j = 2; j < tmppop_size; j++) {
                    if (distance_Matrix[plist->index2][j] < temp1) {
                        temp1 = distance_Matrix[plist->index2][j];
                    }
                }
            } else {
                temp1 = distance_Matrix[plist->index2][0];
                for (j = 1; j < tmppop_size; j++) {
                    if ((j != plist->index) && (distance_Matrix[plist->index2][j] < temp1)) {
                        temp1 = distance_Matrix[plist->index2][j];
                    }
                }
            }
            if (plist->index2 == 0) {
                temp2 = distance_Matrix[plist->index][plist->index2 + 1];
                for (j = 2; j < tmppop_size; j++) {
                    if (distance_Matrix[plist->index][j] < temp2) {
                        temp2 = distance_Matrix[plist->index][j];
                    }
                }
            } else {
                temp2 = distance_Matrix[plist->index2][0];
                for (j = 1; j < tmppop_size; j++) {
                    if ((j != plist->index2) && (distance_Matrix[plist->index][j] < temp2)) {
                        temp2 = distance_Matrix[plist->index][j];
                    }
                }
            }
            temp_list = elist->child;
            if (temp1 < temp2) {
                while (temp_list != NULL) {
                    if (temp_list->index != plist->index2) {
                        temp_list = temp_list->child;
                    } else {
                        for (j = 0; j < tmppop_size; j++) {
                            distance_Matrix[j][plist->index2] = INF;
                        }
                        temp_list = del(temp_list);
                        break;
                    }
                }
            } else {
                while (temp_list != NULL) {
                    if (temp_list->index != plist->index) {
                        temp_list = temp_list->child;
                    } else {
                        for (j = 0; j < tmppop_size; j++) {
                            distance_Matrix[j][plist->index] = INF;
                        }
                        temp_list = del(temp_list);
                        break;
                    }
                }
            }
        }
    }
    temp_list = elist->child;
    i = 0;
    while (temp_list != NULL) {
        copy_ind(&(tmppop->ind[temp_list->index]), &(pop2->ind[i]));
        i++;
        temp_list = temp_list->child;
    }
    while (elist != NULL) {
        temp_list = elist;
        elist = elist->child;
        free(temp_list);
    }
    for (i = 0; i < tmppop_size; i++) {
        free(distance_Matrix[i]);
    }
    free(distance_Matrix);
}

void selection_spea2 (population_real *mixed_pop,int total_size,population_real *archive,int archive_size)
{
    int i, j;
    int num_nondominated;
    individual_real *ind, *ind1, *ind2, *temp_ind;
    population_real *temp_pop;
    int flag;
    double min, temp;

    temp_ind = (individual_real *) malloc (sizeof(individual_real));
    allocate_memory_ind (temp_ind);

    num_nondominated = 0;
    for (i = 0; i < total_size; i++)
    {
        ind = &(mixed_pop->ind[i]);
        if (ind->fitness < 1.0)
        {
            num_nondominated += 1;
        }
    }
    if (num_nondominated <= archive_size)
    {
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
                    min = temp;
                }
            }
            copy_ind (&(mixed_pop->ind[i]), temp_ind);
            copy_ind (&(mixed_pop->ind[flag]), &(mixed_pop->ind[i]));
            copy_ind (temp_ind, &(mixed_pop->ind[flag]));

            copy_ind (&(mixed_pop->ind[i]), &(archive->ind[i]));
        }
    }
    else
    {
        j = 0;
        temp_pop = (population_real *) malloc (sizeof(population_real));
        allocate_memory_pop (temp_pop, num_nondominated);
        for (i = 0; i < total_size; i++)
        {
            ind = &(mixed_pop->ind[i]);
            if (ind->fitness < 1.0)
            {
                copy_ind (&(mixed_pop->ind[i]), &(temp_pop->ind[j]));
                j += 1;
            }
        }
        truncate_pop (temp_pop, num_nondominated, archive, archive_size);
        deallocate_memory_pop (temp_pop, num_nondominated);
        free (temp_pop);
    }
    deallocate_memory_ind (temp_ind);
    free (temp_ind);
}









