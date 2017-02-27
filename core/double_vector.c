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

#include "../header/double_vector.h"
void double_vector_pushback(struct double_vector * head, double value)
{

    print_error(head == NULL,1,"NULL head in double_vector_push\n");
    struct double_vector * ptr = head;
    while(ptr->next != NULL) ptr = ptr-> next;
    struct double_vector *temp;
    temp = (struct double_vector *)malloc(sizeof(struct double_vector));
    ptr -> next = temp;
    temp -> value = value;
    temp -> next = NULL;
}

double double_vector_pop(struct double_vector * head)
{

    print_error(head == NULL,1,"NULL head in double_vector_pop\n");
    struct double_vector * ptr = head;
    while(ptr->next != NULL && ptr->next -> next != NULL) ptr = ptr-> next;
    double value = ptr->next->value;
    free(ptr->next);
    ptr->next = NULL;
    return value;
}

double double_vector_get(struct double_vector * head, int index)
{
    print_error(head == NULL,1,"NULL head in double_vector_get\n");
    struct double_vector * ptr = head;
    int i;
    for(i = 0 ; i < index; i++)
    {
        if(ptr->next != NULL)
            ptr = ptr -> next;
        else return nan("1");
    }
    return ptr ->value;
}

void double_vector_free(struct double_vector * head)
{
    if(head!=NULL)
    {
        double_vector_free(head->next);
    }
    free(head);
}

void double_vector_print(struct double_vector * head)
{
    int index = 0;
    struct double_vector* ptr = head;
    printf("\n");
    while(ptr!=NULL)
    {
        printf("(%d,%lf)",index,ptr->value);
        index ++;
        ptr = ptr->next;
    }
    printf("\n");
    return ;
}