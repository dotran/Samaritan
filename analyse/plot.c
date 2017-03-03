/*
 * analyse.c:
 *  This file contains the functions to record the results, including the population and metric values for analysis.
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
#include "../header/analyse.h"

void plot(char *filename,char *title)
{

    FILE * gnuplot = NULL;
    gnuplot = popen("gnuplot  -persistent","w");

    char command[BUFSIZE_L];

    print_error(gnuplot == NULL, 1, "Error: gnuplot not found.");

    if (number_objective == 2)
    {
        sprintf(command,"set title \"%s\"\n",title);
        fprintf(gnuplot, "%s \n", command);
        sprintf(command,"plot \'%s\'\n",filename);
        fprintf(gnuplot, "%s \n", command);
        fflush(gnuplot);
    }
    else if(number_objective == 3)
    {
        sprintf(command,"set title \"%s\"\n",title);
        fprintf(gnuplot, "%s \n", command);
        fprintf(gnuplot,"set view 60,130 \n");
        fprintf(gnuplot,"set ticslevel 0\n");
        fprintf(gnuplot,"set xrange [0:]\n");
        fprintf(gnuplot,"set yrange [0:]\n");
        fprintf(gnuplot,"set zrange [0:]\n");
        sprintf(command,"splot \'%s\'\n",filename);
        fprintf(gnuplot, "%s \n", command);

        fflush(gnuplot);
    }
    else
    {
        print_error(1,1,"Error: Unsupported number of objective in plot\n");
    }
}