//
// Created by rxc332 on 17-2-22.
//

#include "igd.h"
#include "../header/global.h"
double ** PF_data=NULL;
int number_of_point;
void igd (void *ptr,int id)
{
    int i,j;
    FILE * PF = NULL;
    char PF_name[30];
    double igd_value=0;
    population_real *pop = (population_real *)ptr;
    if (record==NULL)
    {
        record = (double *) malloc ((1+max_generations)*sizeof(double));
        for(i=0;i<1+max_generations;i++)
        {
            record[i]=nan("1");
        }

    }
    if (record_all_run == NULL)
    {
        record_all_run =(double *) malloc ((1+run_index_end)*sizeof(double));
    }
    // calculate IGD

    sprintf(PF_name,"PF/%s_%dD.pf",test_problem,number_objective);

    PF=fopen ( PF_name , "r" );
    if (PF==NULL)
    {
        print_information(EE,2,"Fail to open PF:",PF_name);
    }

    fscanf(PF,"%d",&number_of_point);
    PF_data = (double **) malloc(number_of_point*sizeof(double*));
    for (i=0;i<number_of_point;i++)
        PF_data[i] = (double*) malloc(number_objective*sizeof(double));
    for (i=0;i<number_of_point;i++)
    {
        for (j = 0; j < number_objective; j++)
        {
            fscanf(PF, "%lf", &PF_data[i][j]);
        }
    }

    for(i=0;i<popsize;i++)
    {
        double min_dis = point_distance(pop->ind[i].obj,PF_data[0],number_objective);
        for (j=1;j<number_of_point;j++)
        {
            double cur_dis = point_distance(pop->ind[i].obj,PF_data[j],number_objective);
            if(min_dis>cur_dis)
                min_dis = cur_dis;
        }
        igd_value += min_dis;
    }
    igd_value = igd_value / popsize;
    //
    if (id==max_generations)
        record_all_run[run_index]=igd_value;
    record[id]= igd_value; // dummy

}

void print_global_igd(char *file_name)
{
    int i;
    FILE *fpt;
    fpt = fopen(file_name,"w");
    for(i=run_index_begin;i<=run_index_end;i++)
    {
        fprintf(fpt, "%lf\n",record_all_run[i]);
    }
    free(record_all_run);
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
    record = NULL;
    for(i = 0;i<number_of_point;i++)
        free(PF_data[i]);
    free(PF_data);
    PF_data = NULL;
}

double point_distance(double *a, double *b, int dim)
{
    int i;
    double re = 0;
    for(i=0;i<dim;i++)
    {
        re = re + (a[i]-b[i])*(a[i]-b[i]);
    }
    return sqrt(re);
}