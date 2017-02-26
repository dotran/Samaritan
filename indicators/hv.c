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
# include "./WFG_1.15/wfg.h"
static double *record         = NULL;
static double *record_all_run = NULL;

/* Calculate the IGD value of a population */
void record_hv (void *ptr, int id)
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

    // calculate hv
    record[id] = calculate_hv (ptr);

    if (id == max_generations)
        record_all_run[run_index] = record[id];

    return;
}

/* Calculate HV value */
double calculate_hv (void *ptr)
{
    int i, j;
    double hv_value;
    double min_distance, cur_distance;

    hv_value = hv_wfg(ptr);

    return hv_value;
}

void print_hv (char *file_name)
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

void print_global_hv (char *file_name)
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
