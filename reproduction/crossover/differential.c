/*
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
#include "../../header/reproduction.h"
#include "../../header/global.h"
#include "../../header/rand.h"
void differential (individual_real **parents,individual_real *child)
{
    int i;
    int r;
    double value;
    double yl;

    double yu;

    r = rnd(0, number_variable-1);
    for (i = 0 ; i < number_variable ;i ++)
    {
        yl = variable_lowerbound[i];
        yu = variable_upperbound[i];
        if(rndreal(0,1)<DEFAULT_CR || i == r)
        {
            value = parents[2]->xreal[i] + DEFAULT_F * (parents[0] -> xreal[i] - parents[1]->xreal[i]);
            value = (value>yu)?yu:(value<yl)?yl:value;
        }
        else
        {
            value = parents[2]->xreal[i];
        }
        child->xreal[i] = value;
    }
}