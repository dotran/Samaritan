//
// Created by Renzhi Chen on 2017/2/13.
//

#ifndef CEA_ANALYSE_H
#define CEA_ANALYSE_H

#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../header/global.h"
#include "../header/population.h"
#include "../header/print.h"

#include "igd.h"

void analyse (void *ptr, int id);

int analyse_list[100];
int runtime_output;
int output_interval;

enum analyse_name{VAR,FUN,IGD,END};
static void _mkdir(const char *dir);

void analyse (void *ptr, int id)
{
    int i, j;

    char output_dir_level1[200];
    char output_dir_level2[200];    //lower dir
    char output_file[200];
    char id_char[10];

    sprintf (id_char,"%d",id);
    // set the output dir
    sprintf (output_dir_level1,"./%s_%d_%d/%s/",
            test_problem,
            number_variable,
            number_objective,
            algorithm_name
    );
    sprintf (output_dir_level2,"./%s_%d_%d/%s/%d/",
            test_problem,
            number_variable,
            number_objective,
            algorithm_name,
            number_runs
    );
    _mkdir (output_dir_level2);

    // first time output, init the parameter and output var and fun
    if (id == 1)
    {
        // set analyse list
        for (i = 0; i < 100; i++)
            analyse_list[i] = 0;
        int  read_ptr = 0;
        char name[20];
        while (1)
        {
            int name_c = 0;
            while (analyse_stream[read_ptr] != ' '
                   && analyse_stream[read_ptr] != '\t'
                   && analyse_stream[read_ptr] != '\n'
                   && analyse_stream[read_ptr] != 0)
            {
                name[name_c] = analyse_stream[read_ptr];
                name_c++;
                read_ptr++;
            }
            if (analyse_stream[read_ptr] == 0)
                name[name_c] = 0;
            name[name_c] = 0;

            if (!strcmp(name, "VAR"))
                analyse_list[VAR] = 1;
            else if (!strcmp(name, "FUN"))
                analyse_list[FUN] = 1;
            else if (!strcmp(name, "IGD"))
                analyse_list[IGD] = 1;

            if (analyse_stream[read_ptr] == 0)
                break;
            read_ptr++;
        }
    }

    if (runtime_output == 1 && (id % output_interval == 0|| id == 1 || id == max_generations))
    {
        if (analyse_list[VAR])
        {
            sprintf (output_file, "%s%s.VAR", output_dir_level2, id_char);
            print_variable (output_file, ptr);
        }
        if (analyse_list[FUN])
        {
            sprintf (output_file, "%s%s.FUN", output_dir_level2, id_char);
            print_objective (output_file, ptr);
        }
        if (analyse_list[IGD]) {
            igd (ptr, id);
        }

    }

    if(id == max_generations)
    {
        if (analyse_list[VAR])
        {
            sprintf (output_file, "%s%d.VAR", output_dir_level1, number_runs);
            print_variable (output_file, ptr);
        }
        if (analyse_list[FUN])
        {
            sprintf (output_file, "%s%d.FUN", output_dir_level1, number_runs);
            print_objective (output_file, ptr);
        }
        if (analyse_list[IGD])
        {
            if (runtime_output != 1)
            {
                igd (ptr, id);
                sprintf (output_file, "%s%d.IGD", output_dir_level1, number_runs);
                print_igd (output_file);
            }
        }
    }
}

static void _mkdir (const char *dir)
{
    char tmp[256];
    char *p = NULL;
    size_t len;

    snprintf (tmp, sizeof(tmp), "%s", dir);
    len = strlen (tmp);
    if (tmp[len - 1] == '/')
        tmp[len - 1] = 0;
    for (p = tmp + 1; *p; p++)
    {
        if (*p == '/')
        {
            *p = 0;
            mkdir (tmp, S_IRWXU);
            *p = '/';
        }
    }
    mkdir (tmp, S_IRWXU);
}

#endif // CEA_ANALYSE_H
