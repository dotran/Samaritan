/*
 * initialization.c:
 *  This file contains the functions to perform initialization operations, mostly for reading parameters.
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

# include <string.h>
# include "time.h"
# include "../header/global.h"
# include "../header/rand.h"
# include "../header/print.h"

char algorithm_name[50];
char test_problem[50];
char dummy[50];
char analyse_stream[200];
int run_index;
int run_index_begin;
int run_index_end;
int init_real (char* argv)
{
    int i;
    int random;
    char configFileName[20];

    algorithm_index = 1;

    strcpy (configFileName, argv);
    FILE * config = NULL;
    config = fopen (configFileName, "r");
    if(config == NULL)
    {
        print_information (EE, 2, "Fail to read configure file:", configFileName);
        exit (-1);
    }

    // read from configure file
    fscanf (config, "%s %s", dummy, algorithm_name);
    fscanf (config, "%s %s", dummy, test_problem);
    fscanf (config, "%s %d", dummy, &number_variable);
    fscanf (config, "%s %d", dummy, &number_objective);
    fscanf (config, "%s %d", dummy, &popsize);
    fscanf (config, "%s %d", dummy, &max_generations);
    fscanf (config, "%s %d", dummy, &runtime_output);
    fscanf (config, "%s %d", dummy, &output_interval);
    fscanf (config, "%s %d", dummy, &run_index_begin);
    fscanf (config, "%s %d", dummy, &run_index_end);
    fgets (analyse_stream, 200, config);
    fgets (analyse_stream, 200, config);
    // SBX parameter settings
    pcross_real = 0.9;
    eta_c       = 15.0;

    // polynomial mutation parameter settings
    pmut_real   = 1.0 / number_variable;
    eta_m       = 20.0;

    // initialize a random seed
    srand ((unsigned) time (NULL));
    random = rand () % 1000;
    seed   = (float) random / 1000.0;
    if (seed <= 0.0 || seed >= 1.0)
    {
        printf ("\n Entered seed value is wrong, seed value must be in (0,1) \n");
        exit (1);
    }
    /* DEMO this should read from file*/

    // boundary settings
    variable_lowerbound = (double *)malloc(number_variable * sizeof(double));
    variable_upperbound = (double *)malloc(number_variable * sizeof(double));
    if (!strcmp(test_problem, "ZDT4"))
    {
        variable_lowerbound[0] = 0.0;
        variable_lowerbound[0] = 1.0;
        for (i = 1; i < number_variable; i++)
        {
            variable_lowerbound[i] = -5.0;
            variable_upperbound[i] = 5.0;
        }
    }
    else
    {
        for (i = 0; i < number_variable; i++)
        {
            variable_lowerbound[i] = 0.0;
            variable_upperbound[i] = 1.0;
        }
    }

    return 0;
}
