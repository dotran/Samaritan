//
// Created by rxc332 on 17-3-15.
//
#include "../../header/selection.h"

int* idx ;
double* nicheCount;
struct double_with_index** solMatrix;
double** distMatrix;
double** fitnessMatrix;
struct double_with_index** subpMatrix;
int* statusWoman;
int* next;


int prefers(int x, int y, struct double_with_index* womanPref, int size)
{
    int i;
    for ( i = 0; i < size; i++) {
        int pref = womanPref[i].idx;
        if (pref == x)
            return 1;
        if (pref == y)
            return 0;
    }
    return 0;
}

void stableMatching(int *statusMan,int * statusWoman ,int * next, struct double_with_index** man_pref, struct double_with_index** woman_pref, int menSize, int womenSize)
{

    // Indicates the mating status


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
        int w = man_pref[m][next[m]].idx;
        next[m]++;
        if (statusWoman[w] == NOT_ENGAGED)
        {
            statusMan[m]   = w;
            statusWoman[w] = m;
        }
        else
        {
            int m1 = statusWoman[w];
            if (prefers(m, m1, woman_pref[w], menSize))
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

