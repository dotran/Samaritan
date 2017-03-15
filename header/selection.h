/*
 * selection.h:
 *  This is the header file for the environmental selection operations.
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

# ifndef Samaritan_SELECTION_H
# define Samaritan_SELECTION_H

# include "../header/global.h"
# include "../header/utility.h"
# include "../header/print.h"
# include "../header/rank_sort.h"
# include "../header/population.h"

/* NSGA-II */
void fill_nondominated_sort (population_real* new_pop, population_real* mixed_pop);
void crowding_fill (population_real *mixed_pop, population_real *new_pop, int count, int front_size, list *elite);

void assign_crowding_distance (population_real *pop, int *dist, int **obj_array, int front_size);
void assign_crowding_distance_list (population_real *pop, list *lst, int front_size);
void assign_crowding_distance_indices (population_real *pop, int c1, int c2);

void quicksort_front_obj(population_real *pop, int objcount, int obj_array[], int obj_array_size);
void q_sort_front_obj(population_real *pop, int objcount, int obj_array[], int left, int right);
void quicksort_dist(population_real *pop, int *dist, int front_size);
void q_sort_dist(population_real *pop, int *dist, int left, int right);

/* MOEA/D */
void moead_free ();
void initialize_uniform_weight ();
void initialize_neighborhood ();
void set_weight (double *v, double unit, double sum, int dim);
void update_neighborhood (population_real* pop, individual_real *individual, int subProblemId, int neighborType);
double fitnessFunction (individual_real* individual, double* lambda);
void update_subproblem (population_real* pop, individual_real* individual, int subProblemId, int neighborType);

int* tour_selection (int depth);
void comp_utility (population_real* pop,population_real* saved_values);

/* MOEA/D-STM */

int prefers(int x, int y, struct double_with_index* womanPref, int size);
void stableMatching(int *statusMan,int * statusWoman , int * next, struct double_with_index** man_pref, struct double_with_index** woman_pref, int menSize, int womenSize);
double norm_vector(double* z);
double calculateDistance2(individual_real* individual, double* lambda);

# endif // Samaritan_SELECTION_H
