//
// Created by rxc332 on 17-2-27.
//

#ifndef SAMARITAN_DOUBLE_VECTOR_H
#define SAMARITAN_DOUBLE_VECTOR_H
# include "../header/print.h"
# include "../header/global.h"
# include "../header/double_vector.h"

typedef struct double_vector{
    double value;
    struct  double_vector * next;
}d_vector;

void double_vector_pushback(struct double_vector * head, double value);
double double_vector_pop(struct double_vector * head);
double double_vector_get(struct double_vector * head, int index);
void double_vector_free(struct double_vector * head);
void double_vector_print(struct double_vector * head);
#endif //SAMARITAN_DOUBLE_VECTOR_H
