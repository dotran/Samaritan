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

# include "../header/metaheuristics.h"
/*
 * fillnds.c:
 *  This file contains the functions to perform non-dominated sorting in NSGA-II.
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

population_real *candidate_pop;
population_real* previousExtreme_pop;
population_real* extreme_pop;
double* maxObjValues ;
double* nadirPoint ;
double* u ;
double ** Z ;
int selected_count;
int candidate_count;
int * candidate_flag;
int * association_flag;
double ** distanceMatrix;
double *previousMaxValues;
double * intercepts;
double ** unitDirections;
double * maxArr;
int *associationCount;
int *lastAssociationCount;
int * ReferenceDirection;
double *PerpendicularDistance;
//int * validReferenceDirection;
int maxRank;
int num_candidates;
list ** fronts;
int * fronts_size;
int RemainingCount;
int lastRank;
void associate(population_real *pop, int size)
{

    double d1, d2, lam;
    int i,j,k;

    // Reset Reference Directions Information


    /* Perpendicular distance calculation */
    for (i = 0; i < number_weight; i++)
    {
        for (j = 0; j < size; j++)
        {
            d1 = 0.0;
            lam = 0.0;
            for (k = 0; k < number_objective; k++)
            {
                d1 += (pop->ind[j].obj[k]-ideal_point[k]) * lambda[i][k] / intercepts[k];
                lam += lambda[i][k] * lambda[i][k];
            }
            lam = sqrt(lam);
            d1 = d1 / lam;
            d2 = 0.0;
            for ( k = 0; k < number_objective; k++)
            {
                d2 += pow(
                        ((pop->ind[j].obj[k]-ideal_point[k]) / intercepts[k]
                         - d1 * lambda[i][k]
                         / lam), 2.0);
            }

            // Store the distance in the matrix and in the individual object
            distanceMatrix[j][i] = sqrt(d2);
        }
    }


// Find the closest reference direction to each individual
    for ( i = 0; i < size; i++)
    {
        double minDistance = distanceMatrix[i][0];
        double * refDir = lambda[0];
        int idx = 0;
        for ( j = 1; j < number_weight; j++)
        {
            if (distanceMatrix[i][j] < minDistance)
            {
                idx = j;
                minDistance = distanceMatrix[i][j];
                refDir = lambda[j];
            }
        }
        ReferenceDirection[i] = idx;
        PerpendicularDistance[i] = minDistance;
        //validReferenceDirection[i] = 1;

        //refDir.surroundingIndividuals.add(individuals[i]);
        //individuals[i].validReferenceDirection = true;
    }
    //return distanceMatrix;
}


void assign_rank ( population_real *pop, int size)
{

    //printf("\n========== assign rank ===========\n");
    int i, j;
    int flag;
    int end;
    int front_size = 0;
    int archieve_size = 0;
    int rank = 0;



    candidate_count = 0;

    list *pool;
    list *elite;
    list *temp1, *temp2,*temp3;

    pool  = (list *)malloc(sizeof(list));
    elite = (list *)malloc(sizeof(list));
    pool->index   = -1;
    pool->parent  = NULL;
    pool->child   = NULL;
    elite->index  = -1;
    elite->parent = NULL;
    elite->child  = NULL;

    temp1 = pool;
    for (i = 0; i < size; i++)
    {
        insert (temp1,i);
        temp1 = temp1->child;
    }
    i = 0;
    do
    {
        temp1 = pool->child;
        insert (elite, temp1->index);
        front_size = 1;
        temp2 = elite->child;
        temp1 = del (temp1);
        temp1 = temp1->child;
        do
        {
            temp2 = elite->child;
            if (temp1 == NULL)
                break;

            do
            {
                end  = 0;
                flag = check_dominance (&(pop->ind[temp1->index]), &(pop->ind[temp2->index]));
                if (flag == 1)
                {
                    insert (pool, temp2->index);
                    temp2 = del (temp2);
                    front_size--;
                    temp2 = temp2->child;
                }
                if (flag == 0)
                {
                    temp2 = temp2->child;
                }
                if (flag == -1)
                {
                    end = 1;
                }
            }
            while (end != 1 && temp2 != NULL);
            if (flag == 0 || flag == 1)
            {
                insert (elite, temp1->index);
                front_size++;
                temp1 = del (temp1);
            }
            temp1 = temp1->child;
        }
        while (temp1 != NULL);

        ////////////////////debug
        //printf("\nelite:");
        //temp2 = elite->child;
        //do{
        //    printf("%d ",temp2->index);
        //    temp2=temp2->child;
        //}while(temp2!=NULL);
        //printf("\n");
        ///////////////debug end

        // copy each level into candidate
        if ( candidate_count < size)
        {
            fronts_size[rank] = front_size;

            if(fronts[rank]!=NULL)// need to clean up
            {
                //////// debug
                //printf("rank%d:",rank);
                //temp3 = fronts[rank]->child;
                //do{
                //    printf("(%d)",temp3->index);
                //    temp3 = temp3->child;
                //}while(temp3!=NULL);
                //printf("\n");
                /////// debug end

                temp3 = fronts[rank]->child;
                do
                {
                    //printf("del:%d\n",temp3->index);
                    temp3 = del (temp3);
                    temp3 = temp3->child;
                }
                while (temp3 !=NULL);

                //fronts[rank] = NULL;
            }
            else
            {
                fronts[rank] = malloc(sizeof(list));
            }

            fronts[rank]->index = -1;
            fronts[rank]->parent = NULL;
            fronts[rank]->child = NULL;

            temp3 = fronts[rank];
            maxRank = rank;
            temp2 = elite->child;

            //printf("f[%d]:",rank);
            do
            {
                ////////////////////debug

                //printf("%d ", temp2->index);
                ///////////////debug end

                pop->ind[temp2->index].rank = rank;
                candidate_count += 1;
                insert(temp3,temp2->index);
                temp2 = temp2->child;

            }
            while (temp2 != NULL);
            //printf("(sum:%d)\n",candidate_count);
            rank += 1;
        }

        temp2 = elite->child;
        do
        {
            temp2 = del (temp2);
            temp2 = temp2->child;
        }
        while (elite->child !=NULL);
    }
    while (candidate_count < size);

    // free memory
    while (pool!=NULL)
    {
        temp1 = pool;
        pool = pool->child;
        free (temp1);
    }
    while (elite!=NULL)
    {
        temp1 = elite;
        elite = elite->child;
        free (temp1);
    }
    return;
}
void get_candidate(population_real* mix_pop,population_real *new_pop)
{
    num_candidates = 0;
    int rank = 0;
    list * temp;
    int i;
    while(num_candidates<popsize&&rank<maxRank)
    {
        temp = fronts[rank]->child;
        lastRank = rank;
        selected_count = num_candidates;
        do{

            copy_ind(&mix_pop->ind[temp->index],&candidate_pop->ind[num_candidates]);
            //printf("copy(1) %d(%d) -> %d\n", temp->index, candidate_pop->ind[num_candidates].rank, num_candidates);
            if(num_candidates<popsize) {
                copy_ind(&mix_pop->ind[temp->index], &new_pop->ind[num_candidates]);
                //printf("copy(1) %d(%d) -> %d\n", temp->index, candidate_pop->ind[num_candidates].rank, num_candidates);

            }
            num_candidates = num_candidates +1;
            temp = temp->child;

        }while(temp!=NULL);
        rank  = rank +1;

    }
    //printf("\ncan:%d sel:%d\n",num_candidates,selected_count);

}

/*
void translated_pop(population_real* pop, population_real* nontranslated_pop, int size)
{
    int i,j;
    for(i=0;i<size;i++)
    {
        if(translated[i])
        {
            printf("Translated individuals should NOT be re-translated");
        }
        else
        {
            copy_ind(&pop->ind[i], &nontranslated_pop->ind[i]);
            for(j=0;j<number_objective;j++)
            {
                pop->ind[i].obj[j] = pop->ind[i].obj[j] - ideal_point[j];
            }
        }
    }
}
*/

void getExtremePoints(population_real* pop,int size)
{
    // Calculate unit directions
    int i,j,k;
    for (i = 0; i < number_objective; i++) {
        for (j = 0; j < number_objective; j++) {
            if (i == j) {
                unitDirections[i][j] = 1;
            }
            else {
                unitDirections[i][j] = EPS;     // 1e-6
            }
        }
    }

    //no need to retranslated
    /*
    for ( i = 0; i < number_objective; i++)
    {
            for ( j = 0; j < number_objective; j++) {
                double retranslatedObjValue
                        = previousExtreme_pop->ind[i].obj[j] - (ideal_point[j] - prevIdealPoint[j]);
                previousExtreme_pop->ind[i].obj[j]= retranslatedObjValue;
            }
    }
    */
        // Re-Calculate the previous MAX values of the previous extreme points
        for ( i = 0; i < number_objective; i++) {
            // Set the unit direction (unit direction j)
            double* wDirection = unitDirections[i];
            previousMaxValues[i] = (previousExtreme_pop->ind[i].obj[0]-ideal_point[0]) / wDirection[0];
            for ( k = 1; k < number_objective; k++) {
                double nextValue = (previousExtreme_pop->ind[i].obj[k]-ideal_point[k])/ wDirection[k];
                if (nextValue > previousMaxValues[i]) {
                    previousMaxValues[i] = nextValue;
                }
            }
        }


    for ( i = 0; i < number_objective; i++) {

        double* wDirection = unitDirections[i];


        // Iterate over all the members of the populations
        for ( j = 0; j < size; j++) {
            double max = (pop->ind[j].obj[0]-ideal_point[0]) / wDirection[0];
            for ( k = 1; k < number_objective; k++) {
                double nextValue = (pop->ind[j].obj[k]-ideal_point[k]) / wDirection[k];
                if (nextValue > max) {
                    max = nextValue;
                }
            }
            maxArr[j] = max;
        }
        // Select the minimum value out of maxArr
        int minIndex = 0;
        for ( j = 1; j < size; j++) {
            if (maxArr[j] < maxArr[minIndex]) {
                minIndex = j;
            }
        }

        if ( previousMaxValues[i] < maxArr[minIndex])
        {
            // This means that the previous extreme point was better than
            // the current extreme point and we should retain the previous
            // extreme point instead of replacing it with a new weaker one.
            //extremePoints[i] = previousExtreme_pop[i];
            copy_ind(&previousExtreme_pop->ind[i],&extreme_pop->ind[i]);
        }
        else
        {
            // Now the individual whose index in minIndex in the population is
            // the one representing the extreme factor in the current directions.
            //extremePoints[i] = new Individual(optimizationProblem, individuals[minIndex], individualEvaluator);
            copy_ind(&pop->ind[minIndex],&extreme_pop->ind[i]);
        }
    }
    // Return the extreme points in all basic directions
    //return extremePoints;
}

double* gaussianElimination(double** A, double* b, double *x)
{
    int N = number_objective;
    int p,i,j;
    for ( p = 0; p < N; p++) {

        // find pivot row and swap
        int max = p;
        for (i = p + 1; i < N; i++) {
            if (fabs(A[i][p]) > fabs(A[max][p])) {
                max = i;
            }
        }
        double* temp = A[p];
        A[p] = A[max];
        A[max] = temp;
        double t = b[p];
        b[p] = b[max];
        b[max] = t;

        // singular or nearly singular
        if (fabs(A[p][p]) <= EPS) {
            return NULL;
        }

        // pivot within A and b
        for (i = p + 1; i < N; i++) {
            double alpha = A[i][p] / A[p][p];
            b[i] -= alpha * b[p];
            for ( j = p; j < N; j++) {
                A[i][j] -= alpha * A[p][j];
            }
        }
    }

    // back substitution
    for (i = N - 1; i >= 0; i--) {
        double sum = 0.0;
        for (j = i + 1; j < N; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }
    return x;
}
void niching( population_real * mix_pop,population_real * new_pop)
{
    int i,j;
    list *temp1;
    for(i=0;i<number_weight;i++)
    {
        associationCount[i] = 0;
        association_flag[i] = 1;
        lastAssociationCount[i] = 0;
    }



    for(i=0;i<selected_count;i++)
    {

        associationCount[ReferenceDirection[i]]= associationCount[ReferenceDirection[i]]+1;
    }




    for(i=selected_count;i<num_candidates;i++)
    {

        lastAssociationCount[ReferenceDirection[i]]= lastAssociationCount[ReferenceDirection[i]]+1;
        candidate_flag[i] = 1;
    }



    int remainingIndvsCount = popsize- selected_count;
    struct int_vector* minlist;
    minlist = malloc(sizeof(struct int_vector));
    minlist->next = NULL;
    minlist->value = -1;

    int selected_count2 = 0;
    while(remainingIndvsCount>0)
    {
        /*
        printf("\nremain:%d\n",remainingIndvsCount);
        printf("candidate:");
        for(i=selected_count;i<num_candidates;i++)
            printf("%d",candidate_flag[i]);
        printf("\n");
        for(i=selected_count;i<num_candidates;i++)
            printf("(%d,%lf)",ReferenceDirection[i],PerpendicularDistance[i]);
        printf("\n");
        */
        //printf("\nassociation:");
        //for(i=0;i<number_weight;i++)
        //    printf("%d",association_flag[i]);
        //printf("\nall:");
        //printf("%d",associationCount[0]);
        //for(i=1;i<number_weight;i++)
        //    printf(",%d",associationCount[i]);
        //printf("\nlast");
        //printf("%d",lastAssociationCount[0]);
        //for(i=1;i<number_weight;i++)
        //    printf(",%d",lastAssociationCount[i]);
        //printf("\n");


        int minDirClusterSize = -1;
        for (i = 0; i < number_weight; i++) {
            if(association_flag[i] == 0) continue;
            if (associationCount[i] <= minDirClusterSize||minDirClusterSize==-1) {
                minDirClusterSize = associationCount[i] ;
            }
        }
        int minsize =0;
        //printf("minDirSize:%d, add:",minDirClusterSize);
        minlist = malloc(sizeof(struct int_vector));
        minlist->next = NULL;
        minlist->value = -1;

        for (i = 0; i < number_weight; i++) {
            if(association_flag[i] == 0) continue;
            if (associationCount[i] == minDirClusterSize) {
                int_vector_pushback(minlist,i) ;
                minsize ++;
                //printf("%d ",i);
            }
        }
        //printf("\n");


        int r = rnd(1,minsize);
        int dirIndex = int_vector_get(minlist,r);

        //printf("v:%d(%d) is selected\n",dirIndex,r);
        int_vector_free(minlist);



        if(lastAssociationCount[dirIndex]==0)
        {
            association_flag[dirIndex] = 0;
        }
        else
        {
            int newMemberIndex = -1;

            if(associationCount[dirIndex]==0)
            {
                // select one with smallest
                int minid = -1;
                double minval = -1;
                for(i=selected_count;i<num_candidates;i++)
                {
                    if(dirIndex ==ReferenceDirection[i]&&candidate_flag[i]&&(minid==-1||PerpendicularDistance[i]<minval))
                    {
                        minval = PerpendicularDistance[i];
                        minid = i;
                    }
                }
                newMemberIndex = minid;
            }
            else
            {
                for(i=selected_count;i<num_candidates;i++)
                {
                    if(dirIndex == ReferenceDirection[i]&&candidate_flag[i])
                        break;
                }
                newMemberIndex = i;
            }
            //printf("copy(2)%d: %d(%d) -> %d\n",remainingIndvsCount,newMemberIndex,dirIndex, selected_count+selected_count2);
            copy_ind(&mix_pop->ind[newMemberIndex], &new_pop->ind[selected_count+selected_count2]);

            associationCount[dirIndex] = associationCount[dirIndex]+1;
            lastAssociationCount[dirIndex] = lastAssociationCount[dirIndex] -1;
            candidate_flag[newMemberIndex] = 0;
            selected_count2++;
            remainingIndvsCount --;

        }

    }


}

void getIntercepts(population_real * pop, int size) {
// Calculating the vector of maximum objective values & the Nadir point
// Initialize the structures & set all their initial values to
// Negative Infinity

    int i,j;


    for( i = 0; i< number_objective;i++)
    {
        Z[i] = malloc(number_objective * sizeof(double));
    }

    for ( i = 0; i < number_objective; i++) {
        maxObjValues[i] = -EPS ;
        nadirPoint[i] = -EPS ;
    }
// Traverse all the individuals of the population and get their maximum
// value of objective (The simplest way of calculating the nadir point
// is to get these maximum values among the first front individuals)
    for ( i = 0; i < size; i++)
    {
        for ( j = 0; j < number_objective; j++)
        {
            if (maxObjValues[j] < pop->ind[i].obj[j]-ideal_point[j])
            {
                maxObjValues[j] = pop->ind[i].obj[j]-ideal_point[j];
            }
            if (pop->ind[i].rank == 0)
            {
                if (nadirPoint[j] < pop->ind[i].obj[j]-ideal_point[j])
                {
                    nadirPoint[j] = pop->ind[i].obj[j]-ideal_point[j];
                }
            }
        }
    }


// Caculating the intercepts
 // Create the hyperplane
// Prepare your arrays for gaussian elimination
    for (i = 0; i < number_objective; i++) {
        for ( j = 0; j < number_objective; j++) {
            Z[i][j] = extreme_pop->ind[i].obj[j]-ideal_point[j];
        }
    }

    for (i = 0; i < number_objective; i++) {
        u[i] = 1;
    }
    int useNadir = 0;   // false
// Solve the system of equations using gaussian elimination

    if( gaussianElimination(Z, u,intercepts) == NULL){
    useNadir = 1;
    }

    if (!useNadir) {
        for (i = 0; i < number_objective; i++) {
            intercepts[i] = 1 / intercepts[i];

        }
    }
// If the follwing condition is true this means that you have to resort to the nadir point
    else {
        for (i=0;i<number_objective;i++)
        {
            intercepts[i] = nadirPoint[i];
        }
    }
// If any of the intercepts is still Zero (which means that one of
// the nadir values is Zero), then use the maximum value of each
// objective instead (remember that these values were calculated among
// all the individuals, not just the first-front individuals)
    for (i=0;i<number_objective;i++)
    {
        if (intercepts[i] < EPS)
        {
            for (j = 0; j < number_objective; j++)
                intercepts[j] = maxObjValues[j];
            break;
        }

    }
}


void NSGA3 (population_real* parent_pop, population_real* offspring_pop, population_real* mixed_pop)
{
    int i;
    int generation;

    initialize_uniform_weight ();

    // for assign_rank
    fronts = (list**)malloc(2*popsize*sizeof(list*));
    fronts_size = (int*)malloc(2*popsize*sizeof(int));;
    for(i=0;i<2*popsize;i++)
        fronts[i] = NULL;

    // for gaussian
    maxObjValues = malloc(number_objective * sizeof(double));
    nadirPoint = malloc(number_objective * sizeof(double));
    u = malloc(number_objective * sizeof(double));
    Z = malloc(number_objective * sizeof(double *));

    // for niching
    ReferenceDirection = (int *)malloc(2 * popsize * sizeof(int));
    PerpendicularDistance = (double *)malloc(2 * popsize * sizeof(double));
    association_flag = (int *) malloc(number_weight * sizeof(int));
    candidate_flag = (int *) malloc(2*popsize * sizeof(int));
    intercepts = (double *) malloc( number_objective * sizeof(double));
    previousMaxValues = (double *) malloc( number_objective * sizeof(double));
    maxArr = (double *) malloc( 2*popsize* sizeof(double));
    associationCount = (int *) malloc(number_weight * sizeof(int));
    lastAssociationCount = (int *) malloc(number_weight * sizeof(int));
    distanceMatrix = (double **) malloc( sizeof(double *) * 2 * popsize);
    for(i=0;i<2*popsize;i++)
        distanceMatrix[i] = (double *) malloc(sizeof(double) * popsize);

    unitDirections = (double **) malloc (number_objective * sizeof(double *));
    for(i=0;i<number_objective;i++)
        unitDirections[i] = (double *) malloc(sizeof(double) * popsize);


    previousExtreme_pop = (population_real *) malloc (sizeof(population_real));
    allocate_memory_pop (previousExtreme_pop, number_objective);

    extreme_pop = (population_real *) malloc (sizeof(population_real));
    allocate_memory_pop (extreme_pop, number_objective);

    candidate_pop = (population_real *) malloc (sizeof(population_real));
    allocate_memory_pop (candidate_pop, 2 * popsize);



    generation       = 1;
    evaluation_count = 0;
    printf ("Progress: 1%%");

    initialize_neighborhood ();
    initialize_population_real (parent_pop);
    // population evaluations
    evaluate_population (parent_pop);


    initialize_idealpoint (parent_pop);

    assign_rank(parent_pop,popsize);

    //translated_pop(parent_pop,nontranslated_pop,popsize);

    getExtremePoints(parent_pop, popsize);

    getIntercepts(parent_pop,popsize);

    //resetObjectiveValues(parent_pop, popsize);

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


        assign_rank(mixed_pop,2*popsize);


        for(i=0;i<2*popsize;i++)
            update_ideal_point(&mixed_pop->ind[i]);

        get_candidate(mixed_pop,parent_pop);

        //translated_pop(candidate_pop,nontranslated_pop, num_candidates);

        getExtremePoints(candidate_pop,num_candidates);

        getIntercepts(candidate_pop,num_candidates);

        associate(candidate_pop,num_candidates);

        niching(candidate_pop,parent_pop);
        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, generation, evaluation_count >= max_evaluation);
    }

    // free memory
    free(fronts_size);

    deallocate_memory_pop (candidate_pop, 2 * popsize);
    free (candidate_pop);
    free(associationCount);
    free(lastAssociationCount);
    free(association_flag);
    for(i=0;i<number_objective;i++)
    {
        free(unitDirections[i]);

    }
    for(i=0;i<2*popsize;i++)
        free(distanceMatrix[i]);
    free(candidate_flag);
    free(unitDirections);
    free(distanceMatrix);
    free(fronts);
    free( maxObjValues);
    free(nadirPoint);
    free(u);
    free(Z);

    return;
}
