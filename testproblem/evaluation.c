/*
 * evaluation.c:
 *  This file contains the functions to perform function evaluations.
 *
 * Authors:
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *  Ke Li <k.li@exeter.ac.uk>
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

# include "../header/problems.h"
# include "../header/print.h"

void evaluate_population (population_real* pop)
{
    int i;

    for (i = 0; i < popsize; i++)
        evaluate_individual (&(pop->ind[i]));

    return;
}

void evaluate_individual (individual_real* ind)
{
    int flag;

    flag = 0;
    (strcmp (problem_name, "ZDT1")  != 0)? :(zdt1 (ind->xreal, ind->obj), flag =1);
    (strcmp (problem_name, "ZDT2")  != 0)? :(zdt2 (ind->xreal, ind->obj), flag =1);
    (strcmp (problem_name, "ZDT3")  != 0)? :(zdt3 (ind->xreal, ind->obj), flag =1);
    (strcmp (problem_name, "ZDT4")  != 0)? :(zdt4 (ind->xreal, ind->obj), flag =1);
    (strcmp (problem_name, "ZDT6")  != 0)? :(zdt6 (ind->xreal, ind->obj), flag =1);
    (strcmp (problem_name, "DTLZ1") != 0)? :(dtlz1 (ind->xreal, ind->obj), flag =1);
    (strcmp (problem_name, "DTLZ2") != 0)? :(dtlz2 (ind->xreal, ind->obj), flag =1);
    (strcmp (problem_name, "DTLZ3") != 0)? :(dtlz3 (ind->xreal, ind->obj), flag =1);
    (strcmp (problem_name, "DTLZ4") != 0)? :(dtlz4 (ind->xreal, ind->obj), flag =1);
    (strcmp (problem_name, "DTLZ5") != 0)? :(dtlz5 (ind->xreal, ind->obj), flag =1);
    (strcmp (problem_name, "DTLZ6") != 0)? :(dtlz6 (ind->xreal, ind->obj), flag =1);
    (strcmp (problem_name, "DTLZ7") != 0)? :(dtlz7 (ind->xreal, ind->obj), flag =1);

    print_error (flag == 0, 2, "UNKNOWN test problem: ", problem_name);

    evaluation_count++;

    return;
}
