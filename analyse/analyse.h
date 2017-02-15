//
// Created by Administrator on 2017/2/13.
//

#ifndef CEA_ANALYSE_H
#define CEA_ANALYSE_H


#include "../header/global.h"
#include "../header/population.h"

void analyse(void *ptr, char *fileName);

void analyse(void *ptr, char *fileName)
{
    int i, j;
    FILE *fpt;

    fpt = fopen(fileName,"w");
    population_real *pop;
    pop = ptr;
    for (i = 0; i < popsize; i++)
    {
        /*
        for (j=0; j<x_Dim; j++)
        {
            fprintf(fpt,"%e\t",pop->ind[i].xreal[j]);
        }
        fprintf(fpt,":\t");
        */
        for (j = 0; j < number_objective; j++)
            fprintf (fpt, "%e\t", pop->ind[i].obj[j]);

        fprintf (fpt, "\n");
    }
}

#endif //CEA_ANALYSE_H
