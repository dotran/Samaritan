/*
 * moeadutil.c:
 *  This is the source file for the utility functions for moead.
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

#include "moeadutil.h"
#include "../../header/population.h"
#include "../../header/rank_sort.h"

int C;

void moead_free()
{
    int i;
    free(ideal_point);
    for(i = 0; i< reference_size; i++)
        free(lambda[i]);
    free(lambda);
    for(i = 0 ; i < popsize ; i++)
        free(neighborhood[i]);
    free(neighborhood);
    lambda = NULL;
    ideal_point = NULL;
    neighborhood = NULL;
}


void initialize_uniform_weight () {

    int layer;
    int i, j;
    int * layer_size;
    int * gaps_table;
    int ptr;
    double shrink;
    double * Vec;
    gaps_table = weight_gaps_table[number_objective];
    for (layer = 0; layer < number_objective; layer++)
        if (gaps_table[layer] <= 0)
            break;
    reference_size = 0;
    layer_size = (int *) malloc(sizeof(int) * layer);
    for (i = 0; i < layer; i++) {
        layer_size[i] = combination(number_objective + gaps_table[i] - 1, gaps_table[i]);
        reference_size = reference_size + layer_size[i];
    }

    lambda = (double **) malloc(reference_size * sizeof(double *));
    for (i = 0; i < reference_size; i++)
        lambda[i] = (double *) malloc(number_objective * sizeof(double));
    // free in moead_free

    C = 0;

    ptr = 0;

    shrink = 1;
    int l;
    for (l = 0; l < layer; l++) {
        Vec = (double *) malloc(number_objective * sizeof(double));
        for ( i = 0; i < number_objective; i++)
            Vec[i] = 0;
        setweight( Vec, gaps_table[l], 0, number_objective, number_objective);
        for (i = ptr; i < ptr + layer_size[l]; i++)
            for (j = 0; j < number_objective; j++) {
                lambda[i][j] = lambda[i][j] / gaps_table[l];
                lambda[i][j] = (1 - shrink) / number_objective + shrink * lambda[i][j];
            }
        ptr = ptr + layer_size[l];
        shrink = shrink * 0.8;
        free(Vec);
    }
    free(layer_size);

    FILE * weight_file = fopen("weightvector.out","w");
    for(i=0;i<reference_size;i++)
    {
        for(j =0;j<number_objective;j++)
        {
            fprintf(weight_file,"%lf\t",lambda[i][j]);
        }
        fprintf(weight_file,"\n");
    }
    fclose(weight_file);
}

void setweight( double* v, double unit, double sum, int objdim, int dim)
{

    int i;
    if (dim == objdim) {
        for ( i = 0; i < objdim; i++)
            v[i] = 0;
    }

    if (dim == 1) {

        v[0] = unit - sum;
        for ( i = 0; i < number_objective; i++) {
            fflush(stdout);
            lambda[C][i] = v[i];

        }
        C = C + 1;
        return;
    }
    for ( i = 0; i <= unit - sum; i++) {
        v[dim - 1] = i;
        setweight(v, unit, sum + i, objdim, dim - 1);
    }

}


void initialize_neighborhood() {

    int i,j;
    struct double_s *dis;
    dis = (struct double_s*) malloc(sizeof(struct double_s)*reference_size);

    neighborhood = (int **) malloc(popsize * sizeof(int *));
    // free in moead_free
    for (i = 0; i < popsize; i++)
        neighborhood[i] = (int *) malloc(neighbor_size * sizeof(int));
    for (i = 0; i < popsize; i++) {
        int id = i % reference_size;
        // calculate the distances based on weight vectors
        for (j = 0; j < reference_size; j++) {
            dis[j].x = euclidian_distance(lambda[id], lambda[j],number_objective);
            dis[j].idx = j;
        }
        qsort(dis,reference_size,sizeof(struct double_s),sort_double_cmp);
        for( j = 0 ; j < neighbor_size ; j ++)
        {
            neighborhood[i][j] = dis[j].idx;
        }
    }
    free(dis);
}

void initialize_idealpoint(void * pop)
{
    population_real * ptr = (population_real*) pop;
    int i;
    ideal_point = (double *)malloc(sizeof(double)*number_objective);
    // free in moead_free
    for ( i = 0; i < number_objective; i++) {
        ideal_point[i] = INF;
    }

    for ( i = 0 ;i < popsize ; i ++)
        update_ideal_point(&(ptr->ind[i]));

}


void update_ideal_point(individual_real * p)
{
    int i, n;
    for (n = 0; n < number_objective; n++)
    {
        if (p->obj[n] < ideal_point[n])
        {
            ideal_point[n] = p->obj[n];
        }
    }
}

void update_neighborhood(population_real* pop, individual_real* individual, int subProblemId, int neighborType)
{
    int time,size,i,k;
    double f1, f2;
    int* perm;
    time = 0;
    if (neighborType == NEIGHBOR)
    {
        size = neighbor_size;
    }
    else
    {
        size = reference_size;
    }
    perm= malloc(sizeof(int)*size);

    random_permutation(perm, size);
    for (i = 0; i < size; i++)
    {

        if (neighborType == NEIGHBOR)
        {
            k = neighborhood[subProblemId][perm[i]];
        }
        else
        {
            k = perm[i];
        }


        f1 = fitnessFunction(&(pop->ind[k]), lambda[k]);
        f2 = fitnessFunction(individual, lambda[k]);

        if (f2 < f1 && time < maximumNumberOfReplacedSolutions)
        {
            copy_ind(individual,&(pop->ind[k]));
            time++;
        }
        if (time >= maximumNumberOfReplacedSolutions)
        {
            break;
        }
    }
    free(perm);
    return;
}

double fitnessFunction(individual_real* individual, double* lambda)
{
    double fitness = 0;
    int i;

    if (function_type == TCHE)
    {
        double maxFun = -1.0e+30;

        for (i = 0; i < number_objective; i++)
        {
            double diff = fabs(individual->obj[i] - ideal_point[i]);

            double feval;
            if (lambda[i] == 0)
            {
                feval = 0.0001 * diff;
            }
            else
            {
                feval = diff * lambda[i];
            }
            if (feval > maxFun)
            {
                maxFun = feval;
            }
        }

        fitness = maxFun;
    }
    else if (function_type == AGG)
    {
        double sum = 0.0;
        for (i = 0; i < number_objective; i++)
        {
            sum += (lambda[i]) * individual->obj[i];
        }

        fitness = sum;

    }
    else if (function_type ==PBI)
    {
        double d1, d2, nl;
        double theta = 5.0;

        d1 = d2 = nl = 0.0;

        for (i = 0; i < number_objective; i++)
        {
            d1 += (individual->obj[i] - ideal_point[i]) * lambda[i];
            nl += pow(lambda[i], 2.0);
        }
        nl = sqrt(nl);
        d1 = fabs(d1) / nl;

        for (i = 0; i < number_objective; i++)
        {
            d2 += pow((individual->obj[i] - ideal_point[i]) - d1 * (lambda[i] / nl), 2.0);
        }
        d2 = sqrt(d2);

        fitness = (d1 + theta * d2);
    }
    return fitness;
}