//
// Created by rxc332 on 17-2-21.
//

#ifndef SAMARITAN_IGD_H
#define SAMARITAN_IGD_H
//
// Created by rxc332 on 17-2-20.
//
#include "analyse.h"
#include <float.h>
void igd(void *ptr,int id);
void print_igd(char * file_name);
static double *record = NULL;

void igd(void *ptr,int id)
{
    int i;
    if(record==NULL)
    {
        record = (double *) malloc ((1+max_generations)*sizeof(double));
        for(i=0;i<1+max_generations;i++)
        {
            record[i]=nan("1");
        }

    }
    record[id]= 1.0*id; // dummy

}
void print_igd(char *file_name)
{
    FILE *fpt;
    int i,j;
    fpt = fopen(file_name,"w");
    for(i = 0;i< 1+max_generations;i++) {
        if(!isnan(record[i]))
            fprintf(fpt, "%d\t%lf\n", i,record[i]);
    }
    fclose(fpt);
    free(record);
}
#endif //SAMARITAN_IGD_H
