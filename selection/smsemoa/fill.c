//
// Created by rxc332 on 17-3-22.
//

#include "../../header/global.h"
#include "../../header/dominance.h"
#include "../../header/selection.h"
#include "../../header/population.h"
#include "iwfg.h"

void fill_hv_sort (FILECONTENTS *f,population_real* new_pop, population_real* mixed_pop,int size)
{
    int i, j;
    int flag;
    int end;
    int front_size = 0;
    int archieve_size = 0;
    int rank = 1;

    list *pool;
    list *elite;
    list *temp1, *temp2;

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
                flag = check_dominance (&(mixed_pop->ind[temp1->index]), &(mixed_pop->ind[temp2->index]));
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
        temp2 = elite->child;
        j = i;

        if ( (archieve_size + front_size) <= popsize)
        {
            //printf("\naddF%d:",rank);
            do
            {
                //printf("[%d]->[%d]",temp2->index,i);
                copy_ind (&mixed_pop->ind[temp2->index], &new_pop->ind[i]);
                new_pop->ind[i].rank = rank;
                archieve_size += 1;
                temp2 = temp2->child;
                i += 1;

            }
            while (temp2 != NULL);
            rank += 1;
        }
        else
        {
            //printf("\nselectF%d:\n",rank);
            //temp2 = elite->child;
            //do{
                //printf("[%d]",temp2->index);
            //    temp2 = temp2->child;
            //}while(temp2!=NULL);
            //temp2 = elite->child;
            //do{
            //    for(j=0;j<number_objective;j++)
            //        printf("%lf ",mixed_pop->ind[temp2->index].obj[j]);
            //    printf("\n");
            //    temp2 = temp2->child;
            //}while(temp2!=NULL);
            hv_fill (f,mixed_pop, new_pop, i, front_size, elite);
            archieve_size = popsize;
            for (j = i; j < popsize; j++)
            {
                new_pop->ind[j].rank = rank;
            }
        }
        temp2 = elite->child;
        do
        {
            temp2 = del (temp2);
            temp2 = temp2->child;
        }
        while (elite->child !=NULL);
    }
    while (archieve_size < popsize);

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

void hv_fill (FILECONTENTS *f,population_real *mixed_pop, population_real *new_pop, int count, int front_size, list *elite)
{
    int i, j,k;
    int *dist;
    list *temp;

    dist = (int *)malloc(front_size*sizeof(int));


    i_read_data (f,mixed_pop, elite->child, front_size);

    //i_printContents(f);


    //i_printContents(f);
    i_n = f->fronts[0].n;
    double eh[i_n+2];
    if(i_n == 2)
        i_ihv2(f->fronts[0],eh);
    else
        i_ihv(f->fronts[0],eh);

    //double whole = i_hv(f->fronts[0]);
    //int kk=0;
    //temp = elite->child;
    //do
    //{
    //   i_read_data (f,mixed_pop, elite->child, front_size);
    //    printf("[%d]HVC:%lf/%lf\n",temp->index,whole-i_hv_contribution(f->fronts[0],kk,whole),whole);
    //    kk++;
    //    temp = temp->child;
    //}while(temp!=NULL);
    int id;
    temp = elite->child;
    do{
        id = temp->index;
        int same =0;
        for (j=0;j<number_objective;j++)
        {
            if(fabs((nadir_point[j]-eh[j]) - mixed_pop->ind[id].obj[j]) <1e-4)
                same ++;
        }
        if(same == number_objective)
            break;
        temp = temp->child;


    }while(temp!=NULL);

    //printf("eh[%d]:",id);
    //for( i = 0;i < number_objective ;i ++)
    //{
    //    printf("%lf ",ref_point[i]-eh[i]);
    //}
    //printf("(%lf)\n",eh[number_objective]);
    list * temp2 = elite->child;
    //printf("\nadd[HV]:\n");
    do
    {

        if(temp2->index!=id) {
            //printf("[%d->%d]",temp2->index,count);
            copy_ind(&(mixed_pop->ind[temp2->index]), &(new_pop->ind[count]));
            count++;
        }
        temp2=temp2->child;
    }while(temp2!=NULL);
    free_file_content(f);


    // free memory
    free (dist);

    return;
}
