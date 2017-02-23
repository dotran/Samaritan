/*
 * igd.c:
 *  This file contains the functions to calculate inverted generation distance (IGD).
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

# include "../header/global.h"
# include "../header/indicators.h"
# include "../header/utility.h"

static double *record         = NULL;
static double *record_all_run = NULL;

/* Calculate the IGD value of a population */
void record_igd (void *ptr, int id)
{
    int i;

    if (record == NULL)
    {
        record = (double *) malloc ((max_generations + 1) * sizeof(double));
        for(i = 0; i < max_generations + 1; i++)
            record[i] = nan ("1");
    }
    if (record_all_run == NULL)
        record_all_run = (double *) malloc ((run_index_end + 1) * sizeof(double));

    // calculate IGD
    record[id] = calculate_igd (ptr);

    if (id == max_generations)
        record_all_run[run_index] = record[id];

    return;
}

/* Calculate IGD value */
double calculate_igd (void *ptr)
{
    int i, j;
    double igd_value;
    double min_distance, cur_distance;

    population_real *pop = (population_real *)ptr;

    igd_value = 0.0;
    for (i = 0; i < PF_size; i++)
    {
        min_distance = euclidian_distance (PF_data[i], pop->ind[0].obj, number_objective);
        for (j = 1; j < popsize; j++)
        {
            cur_distance = euclidian_distance (PF_data[i], pop->ind[j].obj, number_objective);
            if(min_distance > cur_distance)
                min_distance = cur_distance;
        }
        igd_value += min_distance;
    }
    igd_value = igd_value / PF_size;

    return igd_value;
}

void print_igd (char *file_name)
{
    int i;
    FILE *fpt;

    fpt = fopen (file_name, "w");
    for (i = 0; i < max_generations + 1; i++)
    {
        if (!isnan(record[i]))
            fprintf (fpt, "%lf\n", record[i]);
    }
    fclose (fpt);

    free (record);
    record = NULL;

    return;
}

void print_global_igd (char *file_name)
{
    int i;
    FILE *fpt;

    fpt = fopen (file_name, "w");
    for (i = run_index_begin; i <= run_index_end; i++)
        fprintf (fpt, "%lf\n", record_all_run[i]);
    fclose (fpt);

    free (record_all_run);
    record_all_run = NULL;

    return;
}
