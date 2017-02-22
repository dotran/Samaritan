//
// Created by Renzhi Chen on 2017/2/13.
//

#ifndef CEA_ANALYSE_H
#define CEA_ANALYSE_H

#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../header/global.h"
#include "../header/population.h"
#include "../header/print.h"

#include "igd.h"

void analyse (void *ptr, int id);

int analyse_list[100];
int runtime_output;
int output_interval;

enum analyse_name{VAR,FUN,IGD,END};

static void _mkdir(const char *dir);
void analyse_all();
void analyse (void *ptr, int id);

#endif // CEA_ANALYSE_H
