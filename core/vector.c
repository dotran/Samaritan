/*
 * moeadutil.c:
 *  This is the source file for the vector struct.
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

#include "../header/vector.h"
/*              int vector                */

int int_vector_size(struct int_vector * head)
{
    struct int_vector * node = head;
    int re =0;
    while(node->next!=NULL)
    {
        re++;
        node = node -> next;
    }
    return re;
}
void int_vector_pushback(struct int_vector * head, int value)
{

    print_error(head == NULL,1,"NULL head in int_vector_push\n");
    struct int_vector * ptr = head;
    while(ptr->next != NULL) ptr = ptr-> next;
    struct int_vector *temp;
    temp = (struct int_vector *)malloc(sizeof(struct int_vector));
    ptr -> next = temp;
    temp -> value = value;
    temp -> next = NULL;
}

int int_vector_pop(struct int_vector * head)
{

    print_error(head == NULL,1,"NULL head in int_vector_pop\n");
    struct int_vector * ptr = head;
    while(ptr->next != NULL && ptr->next -> next != NULL) ptr = ptr-> next;
    int value = ptr->next->value;
    free(ptr->next);
    ptr->next = NULL;
    return value;
}

int int_vector_get(struct int_vector * head, int index)
{
    print_error(head == NULL,1,"NULL head in int_vector_get\n");
    struct int_vector * ptr = head;
    int i;
    for(i = 0 ; i < index; i++)
    {
        if(ptr->next != NULL)
            ptr = ptr -> next;
        else return nan("1");
    }
    return ptr ->value;
}

void int_vector_free(struct int_vector * head)
{
    if(head!=NULL)
    {
        int_vector_free(head->next);
    }
    free(head);
}

void int_vector_print(struct int_vector * head)
{
    int index = 0;
    struct int_vector* ptr = head;
    printf("\n");
    while(ptr!=NULL)
    {
        printf("(%d,%d)",index,ptr->value);
        index ++;
        ptr = ptr->next;
    }
    printf("\n");
    return ;
}

/*              double vector                */

int double_vector_size(struct double_vector * head)
{
    struct double_vector * node = head;
    int re =0;
    while(node->next!=NULL)
    {
        re++;
        node = node -> next;
    }
    return re;
}
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