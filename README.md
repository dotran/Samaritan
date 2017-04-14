# Samaritan
## Change Log

### 2017/03/15

- adding MOEA/D-DRA, MOEA/D-STM and MOEA/D-STM-DRA

### 2017/03/08
- change PLOT from gnuplot to python(mplotlib). Can only plot 2D and 3D objective space.
- add offline plot script plot.py, (example: python plot.py _you_folder_contained_medium_FUN_*.out) 

### 2017/03/03

- combine int_vector and double_vector into vector
- adjust moead structure
- remove openMP option (cannot speedup)
- add PLOT in analyse.c, calling gnuplot(https://sourceforge.net/projects/gnuplot/files/gnuplot/ )

### 2017/02/28
- add MOEAD 
- add openMP option (ATTENTION: Strongly not recomand using openMP with fewer than 8 core)

### 2017/02/27
- change PF format
- add PF for dtlz 1-4 (obj 2-7)
- change from counting generation to counting evaluation
- add DTLZ1-4
- add struct double_vector
