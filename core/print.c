//
// Created by Administrator on 2017/2/13.
//


#include "../header/print.h"
#include "../header/population.h"
#include "../header/global.h"

char symbol[4][3]={"EE","II","DD","dd"};
void printVAR(char *file_name, void * ptr)
{
    FILE *fpt;
    int i,j;
    fpt = fopen(file_name,"w");
    population_real *pop = (population_real*)ptr ;
    for(i = 0;i< popsize;i++) {
        for (j = 0; j < number_variable; j++)
            fprintf(fpt, "%lf\t", pop->ind[i].xreal[j]);
        fprintf(fpt,"\n");
    }
    fclose(fpt);
}
void printFUN(char *file_name, void * ptr)
{
    FILE *fpt;
    int i,j;
    fpt = fopen(file_name,"w");
    population_real *pop = (population_real*)ptr ;
    for(i = 0;i< popsize;i++) {
        for (j = 0; j < number_objective; j++)
            fprintf(fpt, "%lf\t", pop->ind[i].obj[j]);
        fprintf(fpt,"\n");
    }
    fclose(fpt);
}
void printInfo(int level, int n, ...)
{

    int i;
    char * info=NULL;
    va_list vl;
    va_start(vl,n);
    if(DEBUG>=level)
    {
        printf("%s::",symbol[level]);
        for(i=0;i<n;i++)
        {
            info = va_arg(vl,char*);
            printf("%s",info);
        }
        printf("\n");

    }
    va_end(vl);
}


/* print all non-negative args one at a time;
   all args are assumed to be of int type */

