//
// Created by rxc332 on 17-3-10.
//


# include "../header/global.h"
# include "../header/population.h"
# include "../header/reproduction.h"
# include "../header/problems.h"
# include "../header/analyse.h"
# include "../header/initialization.h"
int* idx ;
double* nicheCount;
int**    solPref;
double** solMatrix;
double** distMatrix;
double** fitnessMatrix;
int**    subpPref;
double** subpMatrix;

void QuickSort(double* array, int* idx, int from, int to)
{
    if (from < to) {
        double temp = array[to];
        int tempIdx = idx[to];
        int i = from - 1;
        int j;
        for ( j = from; j < to; j++) {
            if (array[j] <= temp) {
                i++;
                double tempValue = array[j];
                array[j] = array[i];
                array[i] = tempValue;
                int tempIndex = idx[j];
                idx[j] = idx[i];
                idx[i] = tempIndex;
            }
        }
        array[to] = array[i + 1];
        array[i + 1] = temp;
        idx[to] = idx[i + 1];
        idx[i + 1] = tempIdx;
        QuickSort(array, idx, from, i);
        QuickSort(array, idx, i + 1, to);
    }
}


int prefers(int x, int y, int* womanPref, int size)
{
    int i;
    for ( i = 0; i < size; i++) {
        int pref = womanPref[i];
        if (pref == x)
            return 1;
        if (pref == y)
            return 0;
    }
    return 0;
}

void stableMatching(int *statusMan, int** manPref, int** womanPref, int menSize, int womenSize)
{

    // Indicates the mating status

    int* statusWoman = malloc(sizeof(int)*womenSize);
    int* next = malloc(sizeof(int)*womenSize);
    int i;
    int NOT_ENGAGED = -1;

    for (i = 0; i < womenSize; i++)
    {
        statusWoman[i] = NOT_ENGAGED;
        next[i] = 0;
    }

    struct int_vector *freeMen;
    freeMen = malloc(sizeof(struct int_vector));
    freeMen->value = INT_MIN;
    freeMen->next = NULL;

    for (i = 0; i < menSize; i++)
        int_vector_pushback(freeMen,i);

    while (freeMen->next!=NULL)
    {
        int m = int_vector_pop(freeMen);
        int w = manPref[m][next[m]];
        next[m]++;
        if (statusWoman[w] == NOT_ENGAGED)
        {
            statusMan[m]   = w;
            statusWoman[w] = m;
        }
        else
        {
            int m1 = statusWoman[w];
            if (prefers(m, m1, womanPref[w], menSize))
            {
                statusMan[m]   = w;
                statusWoman[w] = m;
                int_vector_pushback(freeMen,m1);
            }
            else
            {
                int_vector_pushback(freeMen,m);
            }
        }
    }

    return;
}


double norm_vector(double* z)
{
    double sum = 0;
    int i;
    for (i = 0; i < number_objective; i++)
        sum += z[i] * z[i];

    return sqrt(sum);
};

double calculateDistance2(individual_real* individual, double* lambda)
{

    int i;
    double distance;
    double distanceSum = 0.0;

    double* vecInd = malloc(sizeof(double)*number_objective);
    double* normalizedObj = malloc(sizeof(double)*number_objective);

    for ( i = 0; i < number_objective; i++)
        distanceSum += individual->obj[i];
    for ( i = 0; i < number_objective; i++)
        normalizedObj[i] = individual->obj[i] / distanceSum;
    for ( i = 0; i < number_objective; i++)
        vecInd[i] = normalizedObj[i] - lambda[i];

    distance = norm_vector(vecInd);

    free(vecInd);
    free(normalizedObj);

    return distance;
}

void stm_selection(population_real* parent_pop,population_real* mixed_pop)
{

    int i,j;

    // Calculate the preference values of solution matrix
    for (i = 0; i < 2*popsize; i++)
    {
        int minIndex = 0;
        for (j = 0; j < popsize; j++) {
            fitnessMatrix[i][j] = fitnessFunction(&(mixed_pop->ind[i]), lambda[j]);
            distMatrix[i][j]  	= calculateDistance2(&(mixed_pop->ind[i]), lambda[j]);
            if (distMatrix[i][j] < distMatrix[i][minIndex])
                minIndex = j;
        }
        nicheCount[minIndex] = nicheCount[minIndex] + 1;
    }

    // calculate the preference values of subproblem matrix and solution matrix
    for (i = 0; i < 2*popsize; i++)
    {
        for (j = 0; j < popsize; j++) {
            subpMatrix[j][i] = fitnessFunction(&(mixed_pop->ind[i]), lambda[j]);
            solMatrix[i][j] = distMatrix[i][j] + nicheCount[j];
        }
    }

    // sort the preference value matrix to get the preference rank matrix
    for ( i = 0; i < popsize; i++)
    {
        for ( j = 0; j < 2*popsize; j++)
            subpPref[i][j] = j;
        QuickSort(subpMatrix[i], subpPref[i], 0, 2*popsize - 1);
    }

    for ( i = 0; i < 2*popsize; i++) {
        for ( j = 0; j < popsize; j++)
            solPref[i][j] = j;
        QuickSort(solMatrix[i], solPref[i], 0, popsize - 1);
    }

    stableMatching(idx,subpPref, solPref, popsize, popsize*2);

    for (i = 0; i < popsize; i++)
    {

        copy_ind(&(mixed_pop->ind[idx[i]]), &(parent_pop->ind[i]));
    }
}


void MOEAD_STM (population_real* parent_pop, population_real* offspring_pop, population_real* mixed_pop)
{
    int i;
    int generation;
    int subproblem_id, neighbor_type;

     idx = malloc(sizeof(int)*popsize);
     nicheCount = malloc(sizeof(double)*popsize);
     solPref   = malloc(sizeof(int*)*popsize*2);
     solMatrix = malloc(sizeof(double*)*popsize*2);
     distMatrix    = malloc(sizeof(int*)*popsize*2);
     fitnessMatrix = malloc(sizeof(int*)*popsize*2);

    for (i = 0; i < 2*popsize; i++)
    {
        solPref[i]   = malloc(sizeof(double)*popsize);
        solMatrix[i] = malloc(sizeof(double)*popsize);
        distMatrix[i]    = malloc(sizeof(double)*popsize);
        fitnessMatrix[i] = malloc(sizeof(double)*popsize);
    }

    subpPref   = malloc(sizeof(int*)*popsize);
    subpMatrix = malloc(sizeof(double*)*popsize);

    for ( i = 0; i < popsize; i++)
    {
        subpPref[i]   = malloc(sizeof(int)*popsize*2);
        subpMatrix[i] = malloc(sizeof(double)*popsize*2);;
    }

    generation       = 1;
    evaluation_count = 0;
    printf ("Progress: 1%%");


    // initialize population
    initialize_population_real (parent_pop);

    // population evaluations
    evaluate_population (parent_pop);

    initialize_uniform_weight ();

    initialize_neighborhood ();

    initialize_idealpoint (parent_pop);

    initialize_nadirpoint(parent_pop);

    permutation = malloc (popsize * sizeof(int));

    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, generation, 0);

    while (evaluation_count < max_evaluation)
    {
        print_progress (generation);

        random_permutation (permutation,popsize);
        for (i = 0; i < popsize; i++)
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
    for (i = 0; i < 2*popsize; i++)
    {
        free(solPref[i]);
        free(solMatrix[i]);
        free(distMatrix[i]);
        free(fitnessMatrix[i]);
    }


    for (i = 0; i < popsize; i++)
    {
        free(subpPref[i]);
        free(subpMatrix[i]);
    }
    return;
}