cmake_minimum_required(VERSION 3.6)
project(Samaritan)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        analyse/analyse.h

        header/dominance.h
        header/global.h
        header/memory.h
        header/population.h
        header/problems.h
        header/rand.h
        header/reproduction.h
        header/selection.h
        header/metaheuristics.h
        header/print.h

        metaheuristics/nsga2.c

        core/initialization.c
        core/dominance.c
        core/list.c
        core/memory.c
        core/population.c
        core/rand.c
        core/print.c

        reproduction/crossover/crossover_nsga2.c
        reproduction/crossover/sbx.c
        reproduction/mating/tournament_selection.c
        reproduction/mutation/mutation_real.c
        reproduction/mutation/polymut.c

        selection/nsga2/crowddist.c
        selection/nsga2/fillnds.c
        selection/nsga2/qsort.c

        testproblem/evaluation.c
        testproblem/DTLZ/DTLZ1.c
        testproblem/ZDT/ZDT1.c

        main.c)

add_executable(Samaritan ${SOURCE_FILES})