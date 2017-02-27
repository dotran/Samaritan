/*
 * problems.h:
 *  This is the header file for benchmark problems.
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

# ifndef Samaritan_PROBLEMS_H
# define Samaritan_PROBLEMS_H

# include "../header/global.h"

void evaluate_population (population_real* pop);
void evaluate_individual (individual_real* ind);
void zdt1 (double *xreal, double *obj);
void dtlz1 (double *xreal, double *obj);
void dtlz2 (double *xreal, double *obj);
void dtlz3 (double *xreal, double *obj);
void dtlz4 (double *xreal, double *obj);

# endif // Samaritan_PROBLEMS_H
