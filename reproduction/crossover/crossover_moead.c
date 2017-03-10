/*
 * crossover_moead.c:
 *  This file contains the functions to perform crossover operations in MOEAD.
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

# include "../../header/reproduction.h"

void crossover_moead_real (population_real *parent_pop, individual_real *offspring, int sub_problem_id, int *neighbor_type)
{
    individual_real **parents = NULL;
    *neighbor_type = choose_neighbor_type ();

    parent_selection (parent_pop, &parents, sub_problem_id, *neighbor_type, 3);
    de (parents, offspring);

    free(parents);

    return;
}

/* Set the neighorhood type used in MOEA/D */
int choose_neighbor_type ()
{
    int neighbor_type ;
    double r;

    r = rndreal (0, 1);
    if (r < neighborhood_selection_probability)
        neighbor_type = NEIGHBOR; //NEIGHBOR;
    else
        neighbor_type = POPULATION; //POPULATION;

    return neighbor_type ;
}