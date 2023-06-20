#ifndef __MEM_H__  // __MEM_H__
#define __MEM_H__

#include "precision.h"
real ****r4dim(int,int,int,int,int,int,int,int);
void del_r4dim(real****,int,int,int,int,int,int,int,int);
real ***r3dim(int,int,int,int,int,int);
void del_r3dim(real***,int,int,int,int,int,int);

float *f1dim(int,int);
void del_f1dim(float*,int,int);
float **f2dim(int,int,int,int);
void del_f2dim(float**,int,int,int,int);
float ***f3dim(int,int,int,int,int,int);
void del_f3dim(float***,int,int,int,int,int,int);

double *d1dim(int,int);
void del_d1dim(double*,int,int);
double **d2dim(int,int,int,int);
void del_d2dim(double**,int,int,int,int);
double ***d3dim(int,int,int,int,int,int);
void del_d3dim(double***,int,int,int,int,int,int);
double ****d4dim(int,int,int,int,int,int,int,int);
void del_d4dim(double****,int,int,int,int,int,int,int,int);
double *****d5dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h);
void del_d5dim(double *****p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h);
real ******r6dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,int x6l,int x6h);
void del_r6dim(real ******p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,int x6l,int x6h);
double ******d6dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,int x6l,int x6h);
void del_d6dim(double ******p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,int x6l,int x6h);
double *******d7dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,int x6l,int x6h,int x7l,int x7h);
void del_d7dim(double *******p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,int x6l
,int x6h,int x7l,int x7h);

int *i1dim(int,int);
void del_i1dim(int*,int,int);
int **i2dim(int,int,int,int);
void del_i2dim(int**,int,int,int,int);
int ***i3dim(int,int,int,int,int,int);
void del_i3dim(int***,int,int,int,int,int,int);
int ***** i5dim(int,int,int,int,int,int,int,int,int,int);
void del_i5dim(int*****,int,int,int,int,int,int,int,int,int,int);
//
int **i2dim(int,int,int,int*);
void del_i2dim(int**,int,int,int,int*);
int ***i3dim(int,int,int,int*,int,int**);
void del_i3dim(int***,int,int,int,int*,int,int**);
//
int ***i3dim(int,int,int,int,int,int*);
void del_i3dim(int***,int,int,int,int,int,int*);

void ***v2dim(int,int,int,int);
void del_v2dim(void***,int,int,int,int);
void ****v3dim(int,int,int,int,int,int);
void del_v3dim(void****,int,int,int,int,int,int);
//
void ****v3dim(int,int,int,int*,int,int**);
void del_v3dim(void****,int,int,int,int*,int,int**);
//
void ***v2dim(int,int,int,int*);
void del_v2dim(void***,int,int,int,int*);
void ****v3dim(int,int,int,int,int,int*);
void del_v3dim(void****,int,int,int,int,int,int*);

#endif             // __MEM_H__
