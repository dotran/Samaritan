/*
 * crossover_nsga2.c:
 *  This file contains the functions to perform crossover operations in NSGA-II.
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

# include "../../header/global.h"
# include "../../header/reproduction.h"
# include "../../header/rand.h"

int choose_neighbor_type();

void crossover_moead_real (population_real *parent_pop, individual_real* offspring, int sub_problem_id, int* neighbor_type)
{

    *neighbor_type = choose_neighbor_type();
    individual_real **parents = NULL;
    parent_selection(parent_pop,&parents,sub_problem_id, *neighbor_type);
    differential(parents,offspring);
    free(parents);

}


int choose_neighbor_type() {
    double r = rndreal(0,1);
    int neighbor_type ;
    if (r < neighborhood_selection_probability)
    {
        neighbor_type = NEIGHBOR; //NEIGHBOR;
    }
    else
    {
        neighbor_type = POPULATION; //POPULATION;
    }
    return neighbor_type ;
}