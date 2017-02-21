//
// Created by rxc332 on 17-2-20.
//

#ifndef SAMARITAN_PRINT_H
#define SAMARITAN_PRINT_H
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#define EE 0
#define II 1
#define DD 2
#define dd 3
int DEBUG;

void printInfo(int level, int n, ...);
void printVAR(char *dir,void* ptr);
void printFUN(char *dir,void* ptr);

#endif //SAMARITAN_PRINT_H
