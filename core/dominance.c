/*
 * dominance.c:
 *  This file contains the functions to perform the dominance check.
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

#include "../header/global.h"

int check_dominance (individual_real *a, individual_real *b)
{
    int i;
    int flag1;
    int flag2;

    flag1 = flag2 = 0;
    for (i = 0; i < number_objective; i++)
    {
        if (a->obj[i] < b->obj[i])
        {
            flag1 = 1;
        }
        else
        {
            if (a->obj[i] > b->obj[i])
            {
                flag2 = 1;
            }
        }
    }
    if (flag1==1 && flag2==0)
    {
        return (1);
    }
    else
    {
        if (flag1==0 && flag2==1)
        {
            return (-1);
        }
        else
        {
            return (0);
        }
    }
}
