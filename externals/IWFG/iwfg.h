#ifndef I_WFG_H_
#define I_WFG_H_

#include "../../header/global.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BEATS(x,y)   (x > y)
#define WORSE(x,y)   (BEATS(y,x) ? (x) : (y))
#define MIN(a,b) (a < b ? (a) : (b))
#define MAX(a,b) (a > b ? (a) : (b))
#define SLICELIMIT 5
typedef double OBJECTIVE;

typedef struct
{
	OBJECTIVE *objectives;
} POINT;

typedef struct
{
	int nPoints;
	int n;
	POINT *points;
} FRONT;

typedef struct
{
        FRONT sprime;   // reduced front 
        int id;         // index in the original list 
        int k;          // next segment to be evaluated 
        double partial; // volume so far 
        int left;       // left child in the heap 
        int right;      // right child in the heap 
} JOB;

typedef struct
{
	int nFronts;
	FRONT *fronts;
} FILECONTENTS;


typedef struct
{
    double width;
    FRONT front;
    int index;
} SLICE;

void i_printContents(FILECONTENTS *f);
void i_read_data(FILECONTENTS *f,void *ptr,list * l,int size);
void free_file_content(FILECONTENTS *f);

double i_hv_contribution(FRONT ps, int id,double whole);
double i_hv(FRONT ps);
int i_slicingDepth(int d);
void i_ihv2(FRONT ps, double *min);
void i_ihv(FRONT ps, double *min);

int i_n;     // the number of objectives
POINT i_ref; // the reference point
POINT i_dirs;// records the directions of the objectives

FRONT *i_fs;    // memory management stuff

int i_maxm; // maxmimum number of points
int i_maxn; // maximum number of objectives
int i_safe;


double* partial; //partial exclusive hypervolumes
int* heap; //heap-based priority queue
int heapsize; //number of points in queue
SLICE **stacks; //set of slices per point per slicing depth
int *stacksize; //current slicing depth per point

int* gorder; //objective order used by comparison functions
int** torder; //order of objectives per point
int** tcompare;
FRONT* fsorted; //front sorted in each objective


#endif
