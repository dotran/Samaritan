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

#define EE 0
#define II 1
#define DD 2
#define dd 3

# include "time.h"
# include "../header/global.h"
# include "../header/rand.h"
# include "../header/memory.h"
# include "print.h"

int init_real ()
{
    int i;
    int random;

    algorithm_index = 1;
    problem_index   = 1;

//    if(argc == 1)
//    {
//        print(II, 1, "Reading Default Configure File configure.cfg");
//    }
//    else if(argc > 1)
//    {
//        print(II, 2, "Reading from file:",argv[1]);

//    }
//    else
//    {
//        print(II, 2, "Extra args ignored, reading from file:",argv[1]);

//    }

    // initialize a random seed
    srand ((unsigned) time (NULL));
    random = rand () % 1000;
    seed = (float) random / 1000.0;
    if (seed <= 0.0 || seed >= 1.0)
    {
        printf ("\n Entered seed value is wrong, seed value must be in (0,1) \n");
        exit (1);
    }
    /* DEMO this should read from file*/
    popsize          = 100;
    number_variable  = 7;
    number_objective = 3;
    max_generations  = 1000;

    // boundary settings
    variable_lowerbound = (double *)malloc(number_variable * sizeof(double));
    variable_upperbound = (double *)malloc(number_variable * sizeof(double));
    for (i = 0; i < number_variable; i++)
    {
        variable_lowerbound[i] = 0.0;
        variable_upperbound[i] = 1.0;
    }

    // parameter initialize
    pcross_real = 0.9;
    pmut_real   = 1.0 / number_variable;
    eta_c       = 15.0;
    eta_m       = 20.0;

    return 0;
}
