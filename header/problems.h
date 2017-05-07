/*
 * problems.h:
 *  This is the header file for benchmark problems.
 *
 * Authors:
 *  Ke Li <k.li@exeter.ac.uk>
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
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
void zdt2 (double *xreal, double *obj);
void zdt3 (double *xreal, double *obj);
void zdt4 (double *xreal, double *obj);
void zdt6 (double *xreal, double *obj);
void dtlz1 (double *xreal, double *obj);
void dtlz2 (double *xreal, double *obj);
void dtlz3 (double *xreal, double *obj);
void dtlz4 (double *xreal, double *obj);
void dtlz5 (double *xreal, double *obj);
void dtlz6 (double *xreal, double *obj);
void dtlz7 (double *xreal, double *obj);
void uf1 (double *xreal, double *obj);
void uf2 (double *xreal, double *obj);
void uf3 (double *xreal, double *obj);
void uf4 (double *xreal, double *obj);
void uf5 (double *xreal, double *obj);
void uf6 (double *xreal, double *obj);
void uf7 (double *xreal, double *obj);
void uf8 (double *xreal, double *obj);
void uf9 (double *xreal, double *obj);
void uf10 (double *xreal, double *obj);
void c1dtlz1 (double *xreal, double *obj,double *cv);
void c2dtlz2 (double *xreal, double *obj,double *cv);
# endif // Samaritan_PROBLEMS_H
