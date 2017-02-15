//
// Created by Administrator on 2017/2/13.
//

#ifndef CEA_PRINT_H
#define CEA_PRINT_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#define EE 0
#define II 1
#define DD 2
#define dd 3
//extern int DEBUG;

char symbol[4][3]={"EE","II","DD","dd"};

void print(int level, int n, ...)
{

    int i;
    char * info=NULL;
    va_list vl;
    va_start(vl,n);
//    if(DEBUG>=level)
//    {
//        printf("%s::",symbol[level]);
//        for(i=0;i<n;i++)
//        {
//            info = va_arg(vl,char*);
//            printf("%s",info);
//        }
//        printf("\n");
//
//    }
    va_end(vl);
}

#include <stdarg.h>

/* print all non-negative args one at a time;
   all args are assumed to be of int type */


#endif //CEA_PRINT_H
