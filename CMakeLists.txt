cmake_minimum_required(VERSION 3.6)
project(Samaritan)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        header/dominance.h
        header/global.h
        header/memory.h
        header/indicators.h
        header/metaheuristics.h
        header/population.h
        header/print.h
        header/problems.h
        header/rand.h
        header/reproduction.h
        header/selection.h
        header/utility.h
        header/analyse.h
        header/double_vector.h

        core/dominance.c
        core/initialization.c
        core/list.c
        core/memory.c
        core/population.c
        core/print.c
        core/rand.c
        core/utility.c
        core/double_vector.c

        metaheuristics/nsga2.c
        metaheuristics/moead.c

        metaheuristics/util/moeadutil.h
        metaheuristics/util/moeadutil.c

        selection/nsga2/crowddist.c
        selection/nsga2/fillnds.c
        selection/nsga2/qsort.c

        reproduction/crossover/crossover_nsga2.c
        reproduction/crossover/crossover_moead.c
        reproduction/crossover/sbx.c
        reproduction/crossover/differential.c
        reproduction/mating/tournament_selection.c
        reproduction/mating/moead_selection.c
        reproduction/mutation/mutation_real.c
        reproduction/mutation/polymut.c

        analyse/analyse.c

        indicators/hv.c
        indicators/igd.c

        indicators/WFG_1.15/read.c
        indicators/WFG_1.15/wfg.c
        indicators/WFG_1.15/wfg.h


        testproblem/DTLZ/DTLZ1.c
        testproblem/DTLZ/DTLZ2.c
        testproblem/DTLZ/DTLZ3.c
        testproblem/DTLZ/DTLZ4.c
        testproblem/ZDT/ZDT1.c
        testproblem/evaluation.c

        config.txt
        main.c
        )

add_executable(Samaritan ${SOURCE_FILES})

target_link_libraries(Samaritan m)
