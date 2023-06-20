#include "mem.h"

real ****r4dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1,nx4=x4h-x4l+1;
  real ****p;
  p=new real*** [nx1] - x1l;
  p[x1l]=new real** [nx1*nx2] - x2l;
  p[x1l][x2l]=new real* [nx1*nx2*nx3] - x3l;
  p[x1l][x2l][x3l]=new real [nx1*nx2*nx3*nx4] - x4l;
  for(int x3=x3l+1;x3<=x3h;++x3) p[x1l][x2l][x3]=p[x1l][x2l][x3-1]+nx4;
  for(int x2=x2l+1;x2<=x2h;++x2){
    p[x1l][x2]=p[x1l][x2-1]+nx3;
    p[x1l][x2][x3l]=p[x1l][x2-1][x3l]+nx3*nx4;
    for(int x3=x3l+1;x3<=x3h;++x3) p[x1l][x2][x3]=p[x1l][x2][x3-1]+nx4;
  }
  for(int x1=x1l+1;x1<=x1h;++x1) {
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    p[x1][x2l][x3l]=p[x1-1][x2l][x3l]+nx2*nx3*nx4;
    for(int x3=x3l+1;x3<=x3h;++x3) p[x1][x2l][x3]=p[x1][x2l][x3-1]+nx4;
    for(int x2=x2l+1;x2<=x2h;++x2){
      p[x1][x2]=p[x1][x2-1]+nx3;
      p[x1][x2][x3l]=p[x1][x2-1][x3l]+nx3*nx4;
      for(int x3=x3l+1;x3<=x3h;++x3) p[x1][x2][x3]=p[x1][x2][x3-1]+nx4;
    }
  }
  return p;
}

void del_r4dim(real ****p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h)
{
  delete[] (p[x1l][x2l][x3l]+x4l);
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

real ***r3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  real ***p;
  p=new real** [nx1] - x1l;
  p[x1l]=new real* [nx1*nx2] - x2l;
  p[x1l][x2l]=new real [nx1*nx2*nx3] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

void del_r3dim(real ***p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

float *f1dim(int x1l,int x1h)
{
  return new float [x1h-x1l+1] - x1l;
}

void del_f1dim(float *p,int x1l,int x1h)
{
  delete[] (p+x1l);
}
float **f2dim(int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  float **p;
  p=new float* [nx1] - x1l;
  p[x1l]=new float [nx1*nx2] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

void del_f2dim(float **p,int x1l, int x1h, int x2l, int x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

float ***f3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  float ***p;
  p=new float** [nx1] - x1l;
  p[x1l]=new float* [nx1*nx2] - x2l;
  p[x1l][x2l]=new float [nx1*nx2*nx3] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

void del_f3dim(float ***p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

double *d1dim(int x1l,int x1h)
{
  return new double [x1h-x1l+1] - x1l;
}

void del_d1dim(double *p,int x1l,int x1h)
{
  delete[] (p+x1l);
}

double **d2dim(int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  double **p;
  p=new double* [nx1] - x1l;
  p[x1l]=new double [nx1*nx2] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

void del_d2dim(double **p,int x1l, int x1h, int x2l, int x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

double ***d3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  double ***p;
  p=new double** [nx1] - x1l;
  p[x1l]=new double* [nx1*nx2] - x2l;
  p[x1l][x2l]=new double [nx1*nx2*nx3] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

void del_d3dim(double ***p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

double ****d4dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1,nx4=x4h-x4l+1;
  double ****p;
  p=new double*** [nx1] - x1l;
  p[x1l]=new double** [nx1*nx2] - x2l;
  p[x1l][x2l]=new double* [nx1*nx2*nx3] - x3l;
  p[x1l][x2l][x3l]=new double [nx1*nx2*nx3*nx4] - x4l;
  for(int x3=x3l+1;x3<=x3h;++x3) p[x1l][x2l][x3]=p[x1l][x2l][x3-1]+nx4;
  for(int x2=x2l+1;x2<=x2h;++x2){
    p[x1l][x2]=p[x1l][x2-1]+nx3;
    p[x1l][x2][x3l]=p[x1l][x2-1][x3l]+nx3*nx4;
    for(int x3=x3l+1;x3<=x3h;++x3) p[x1l][x2][x3]=p[x1l][x2][x3-1]+nx4;
  }
  for(int x1=x1l+1;x1<=x1h;++x1) {
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    p[x1][x2l][x3l]=p[x1-1][x2l][x3l]+nx2*nx3*nx4;
    for(int x3=x3l+1;x3<=x3h;++x3) p[x1][x2l][x3]=p[x1][x2l][x3-1]+nx4;
    for(int x2=x2l+1;x2<=x2h;++x2){
      p[x1][x2]=p[x1][x2-1]+nx3;
      p[x1][x2][x3l]=p[x1][x2-1][x3l]+nx3*nx4;
      for(int x3=x3l+1;x3<=x3h;++x3) p[x1][x2][x3]=p[x1][x2][x3-1]+nx4;
    }
  }
  return p;
}

void del_d4dim(double ****p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h)
{
  delete[] (p[x1l][x2l][x3l]+x4l);
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

double *****d5dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1,nx4=x4h-x4l+1,nx5=x5h-x5l+1;
  double *****p;
  p=new double**** [nx1] - x1l;
  p[x1l]=new double*** [nx1*nx2] - x2l;
  p[x1l][x2l]=new double** [nx1*nx2*nx3] - x3l;
  p[x1l][x2l][x3l]=new double* [nx1*nx2*nx3*nx4] - x4l;
  p[x1l][x2l][x3l][x4l]=new double [nx1*nx2*nx3*nx4*nx5] - x5l;
  for(int x4=x4l+1;x4<=x4h;++x4) p[x1l][x2l][x3l][x4]=p[x1l][x2l][x3l][x4-1]+nx5;
  for(int x3=x3l+1;x3<=x3h;++x3){
    p[x1l][x2l][x3]=p[x1l][x2l][x3-1]+nx4;    
    p[x1l][x2l][x3][x4l]=p[x1l][x2l][x3-1][x4l]+nx4*nx5;
    for(int x4=x4l+1;x4<=x4h;++x4) p[x1l][x2l][x3][x4]=p[x1l][x2l][x3][x4-1]+nx5;
  }
  for(int x2=x2l+1;x2<=x2h;++x2){
    p[x1l][x2]=p[x1l][x2-1]+nx3;
    p[x1l][x2][x3l]=p[x1l][x2-1][x3l]+nx3*nx4;
    p[x1l][x2][x3l][x4l]=p[x1l][x2-1][x3l][x4l]+nx3*nx4*nx5;
    for(int x4=x4l+1;x4<=x4h;++x4) p[x1l][x2][x3l][x4]=p[x1l][x2][x3l][x4-1]+nx5;
    for(int x3=x3l+1;x3<=x3h;++x3){
      p[x1l][x2][x3]=p[x1l][x2][x3-1]+nx4;
      p[x1l][x2][x3][x4l]=p[x1l][x2][x3-1][x4l]+nx4*nx5;
      for(int x4=x4l+1;x4<=x4h;++x4) p[x1l][x2][x3][x4]=p[x1l][x2][x3][x4-1]+nx5;
    }
  }	
  for(int x1=x1l+1;x1<=x1h;++x1) {
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    p[x1][x2l][x3l]=p[x1-1][x2l][x3l]+nx2*nx3*nx4;
    p[x1][x2l][x3l][x4l]=p[x1-1][x2l][x3l][x4l]+nx2*nx3*nx4*nx5;
    for(int x4=x4l+1;x4<=x4h;++x4) p[x1][x2l][x3l][x4]=p[x1][x2l][x3l][x4-1]+nx5;
    for(int x3=x3l+1;x3<=x3h;++x3){
      p[x1][x2l][x3]=p[x1][x2l][x3-1]+nx4;
      p[x1][x2l][x3][x4l]=p[x1][x2l][x3-1][x4l]+nx4*nx5;
      for(int x4=x4l+1;x4<=x4h;++x4) p[x1][x2l][x3][x4]=p[x1][x2l][x3][x4-1]+nx5;
    }
    for(int x2=x2l+1;x2<=x2h;++x2){
      p[x1][x2]=p[x1][x2-1]+nx3;
      p[x1][x2][x3l]=p[x1][x2-1][x3l]+nx3*nx4;
      p[x1][x2][x3l][x4l]=p[x1][x2-1][x3l][x4l]+nx3*nx4*nx5;
      for(int x4=x4l+1;x4<=x4h;++x4) p[x1][x2][x3l][x4]=p[x1][x2][x3l][x4-1]+nx5;
      for(int x3=x3l+1;x3<=x3h;++x3){
        p[x1][x2][x3]=p[x1][x2][x3-1]+nx4;
        p[x1][x2][x3][x4l]=p[x1][x2][x3-1][x4l]+nx4*nx5;
        for(int x4=x4l+1;x4<=x4h;++x4) p[x1][x2][x3][x4]=p[x1][x2][x3][x4-1]+nx5;
      }
    }
  }
  return p;
}

void del_d5dim(double *****p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h)
{
  delete[] (p[x1l][x2l][x3l][x4l]+x5l);
  delete[] (p[x1l][x2l][x3l]+x4l);
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}
real ******r6dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,int x6l,int x6h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1,nx4=x4h-x4l+1,nx5=x5h-x5l+1,nx6=x6h-x6l+1;
  real ******p;
  p=new real***** [nx1] - x1l;
  p[x1l]=new real**** [nx1*nx2] - x2l;
  p[x1l][x2l]=new real*** [nx1*nx2*nx3] - x3l;
  p[x1l][x2l][x3l]=new real** [nx1*nx2*nx3*nx4] - x4l;
  p[x1l][x2l][x3l][x4l]=new real* [nx1*nx2*nx3*nx4*nx5] - x5l;
  p[x1l][x2l][x3l][x4l][x5l]=new real [nx1*nx2*nx3*nx4*nx5*nx6] - x6l;
  for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2l][x3l][x4l][x5]=p[x1l][x2l][x3l][x4l][x5-1]+nx6;
  for(int x4=x4l+1;x4<=x4h;++x4){
    p[x1l][x2l][x3l][x4]=p[x1l][x2l][x3l][x4-1]+nx5;    
    p[x1l][x2l][x3l][x4][x5l]=p[x1l][x2l][x3l][x4-1][x5l]+nx5*nx6;
    for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2l][x3l][x4][x5]=p[x1l][x2l][x3l][x4][x5-1]+nx6;
  }
  for(int x3=x3l+1;x3<=x3h;++x3){
    p[x1l][x2l][x3]=p[x1l][x2l][x3-1]+nx4;
    p[x1l][x2l][x3][x4l]=p[x1l][x2l][x3-1][x4l]+nx4*nx5;
    p[x1l][x2l][x3][x4l][x5l]=p[x1l][x2l][x3-1][x4l][x5l]+nx4*nx5*nx6;
    for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2l][x3][x4l][x5]=p[x1l][x2l][x3][x4l][x5-1]+nx6;
    for(int x4=x4l+1;x4<=x4h;++x4){
      p[x1l][x2l][x3][x4]=p[x1l][x2l][x3][x4-1]+nx5;    
      p[x1l][x2l][x3][x4][x5l]=p[x1l][x2l][x3][x4-1][x5l]+nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2l][x3][x4][x5]=p[x1l][x2l][x3][x4][x5-1]+nx6;
    }
  }
  for(int x2=x2l+1;x2<=x2h;++x2){
    p[x1l][x2]=p[x1l][x2-1]+nx3;
    p[x1l][x2][x3l]=p[x1l][x2-1][x3l]+nx3*nx4;
    p[x1l][x2][x3l][x4l]=p[x1l][x2-1][x3l][x4l]+nx3*nx4*nx5;
    p[x1l][x2][x3l][x4l][x5l]=p[x1l][x2-1][x3l][x4l][x5l]+nx3*nx4*nx5*nx6;
    for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2][x3l][x4l][x5]=p[x1l][x2][x3l][x4l][x5-1]+nx6;
    for(int x4=x4l+1;x4<=x4h;++x4){
      p[x1l][x2][x3l][x4]=p[x1l][x2][x3l][x4-1]+nx5;    
      p[x1l][x2][x3l][x4][x5l]=p[x1l][x2][x3l][x4-1][x5l]+nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2][x3l][x4][x5]=p[x1l][x2][x3l][x4][x5-1]+nx6;
    }
    for(int x3=x3l+1;x3<=x3h;++x3){
      p[x1l][x2][x3]=p[x1l][x2][x3-1]+nx4;
      p[x1l][x2][x3][x4l]=p[x1l][x2][x3-1][x4l]+nx4*nx5;
      p[x1l][x2][x3][x4l][x5l]=p[x1l][x2][x3-1][x4l][x5l]+nx4*nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2][x3][x4l][x5]=p[x1l][x2][x3][x4l][x5-1]+nx6;
      for(int x4=x4l+1;x4<=x4h;++x4){
        p[x1l][x2][x3][x4]=p[x1l][x2][x3][x4-1]+nx5;    
        p[x1l][x2][x3][x4][x5l]=p[x1l][x2][x3][x4-1][x5l]+nx5*nx6;
        for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2][x3][x4][x5]=p[x1l][x2][x3][x4][x5-1]+nx6;
      }
    }
  }
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    p[x1][x2l][x3l]=p[x1-1][x2l][x3l]+nx2*nx3*nx4;
    p[x1][x2l][x3l][x4l]=p[x1-1][x2l][x3l][x4l]+nx2*nx3*nx4*nx5;
    p[x1][x2l][x3l][x4l][x5l]=p[x1-1][x2l][x3l][x4l][x5l]+nx2*nx3*nx4*nx5*nx6;
    for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2l][x3l][x4l][x5]=p[x1][x2l][x3l][x4l][x5-1]+nx6;
    for(int x4=x4l+1;x4<=x4h;++x4){
      p[x1][x2l][x3l][x4]=p[x1][x2l][x3l][x4-1]+nx5;    
      p[x1][x2l][x3l][x4][x5l]=p[x1][x2l][x3l][x4-1][x5l]+nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2l][x3l][x4][x5]=p[x1][x2l][x3l][x4][x5-1]+nx6;
    }
    for(int x3=x3l+1;x3<=x3h;++x3){
      p[x1][x2l][x3]=p[x1][x2l][x3-1]+nx4;
      p[x1][x2l][x3][x4l]=p[x1][x2l][x3-1][x4l]+nx4*nx5;
      p[x1][x2l][x3][x4l][x5l]=p[x1][x2l][x3-1][x4l][x5l]+nx4*nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2l][x3][x4l][x5]=p[x1][x2l][x3][x4l][x5-1]+nx6;
      for(int x4=x4l+1;x4<=x4h;++x4){
        p[x1][x2l][x3][x4]=p[x1][x2l][x3][x4-1]+nx5;    
        p[x1][x2l][x3][x4][x5l]=p[x1][x2l][x3][x4-1][x5l]+nx5*nx6;
        for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2l][x3][x4][x5]=p[x1][x2l][x3][x4][x5-1]+nx6;
      }
    }
    for(int x2=x2l+1;x2<=x2h;++x2){
      p[x1][x2]=p[x1][x2-1]+nx3;
      p[x1][x2][x3l]=p[x1][x2-1][x3l]+nx3*nx4;
      p[x1][x2][x3l][x4l]=p[x1][x2-1][x3l][x4l]+nx3*nx4*nx5;
      p[x1][x2][x3l][x4l][x5l]=p[x1][x2-1][x3l][x4l][x5l]+nx3*nx4*nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2][x3l][x4l][x5]=p[x1][x2][x3l][x4l][x5-1]+nx6;
      for(int x4=x4l+1;x4<=x4h;++x4){
        p[x1][x2][x3l][x4]=p[x1][x2][x3l][x4-1]+nx5;    
        p[x1][x2][x3l][x4][x5l]=p[x1][x2][x3l][x4-1][x5l]+nx5*nx6;
        for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2][x3l][x4][x5]=p[x1][x2][x3l][x4][x5-1]+nx6;
      }
      for(int x3=x3l+1;x3<=x3h;++x3){
        p[x1][x2][x3]=p[x1][x2][x3-1]+nx4;
        p[x1][x2][x3][x4l]=p[x1][x2][x3-1][x4l]+nx4*nx5;
        p[x1][x2][x3][x4l][x5l]=p[x1][x2][x3-1][x4l][x5l]+nx4*nx5*nx6;
        for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2][x3][x4l][x5]=p[x1][x2][x3][x4l][x5-1]+nx6;
        for(int x4=x4l+1;x4<=x4h;++x4){
          p[x1][x2][x3][x4]=p[x1][x2][x3][x4-1]+nx5;    
          p[x1][x2][x3][x4][x5l]=p[x1][x2][x3][x4-1][x5l]+nx5*nx6;
          for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2][x3][x4][x5]=p[x1][x2][x3][x4][x5-1]+nx6;
        }
      }
    }
  }
  return p;
}

void del_r6dim(real ******p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,int x6l,int x6h)
{
  delete[] (p[x1l][x2l][x3l][x4l][x5l]+x6l);
  delete[] (p[x1l][x2l][x3l][x4l]+x5l);
  delete[] (p[x1l][x2l][x3l]+x4l);
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

double ******d6dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,int x6l,int x6h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1,nx4=x4h-x4l+1,nx5=x5h-x5l+1,nx6=x6h-x6l+1;
  double ******p;
  p=new double***** [nx1] - x1l;
  p[x1l]=new double**** [nx1*nx2] - x2l;
  p[x1l][x2l]=new double*** [nx1*nx2*nx3] - x3l;
  p[x1l][x2l][x3l]=new double** [nx1*nx2*nx3*nx4] - x4l;
  p[x1l][x2l][x3l][x4l]=new double* [nx1*nx2*nx3*nx4*nx5] - x5l;
  p[x1l][x2l][x3l][x4l][x5l]=new double [nx1*nx2*nx3*nx4*nx5*nx6] - x6l;
  for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2l][x3l][x4l][x5]=p[x1l][x2l][x3l][x4l][x5-1]+nx6;
  for(int x4=x4l+1;x4<=x4h;++x4){
    p[x1l][x2l][x3l][x4]=p[x1l][x2l][x3l][x4-1]+nx5;    
    p[x1l][x2l][x3l][x4][x5l]=p[x1l][x2l][x3l][x4-1][x5l]+nx5*nx6;
    for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2l][x3l][x4][x5]=p[x1l][x2l][x3l][x4][x5-1]+nx6;
  }
  for(int x3=x3l+1;x3<=x3h;++x3){
    p[x1l][x2l][x3]=p[x1l][x2l][x3-1]+nx4;
    p[x1l][x2l][x3][x4l]=p[x1l][x2l][x3-1][x4l]+nx4*nx5;
    p[x1l][x2l][x3][x4l][x5l]=p[x1l][x2l][x3-1][x4l][x5l]+nx4*nx5*nx6;
    for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2l][x3][x4l][x5]=p[x1l][x2l][x3][x4l][x5-1]+nx6;
    for(int x4=x4l+1;x4<=x4h;++x4){
      p[x1l][x2l][x3][x4]=p[x1l][x2l][x3][x4-1]+nx5;    
      p[x1l][x2l][x3][x4][x5l]=p[x1l][x2l][x3][x4-1][x5l]+nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2l][x3][x4][x5]=p[x1l][x2l][x3][x4][x5-1]+nx6;
    }
  }
  for(int x2=x2l+1;x2<=x2h;++x2){
    p[x1l][x2]=p[x1l][x2-1]+nx3;
    p[x1l][x2][x3l]=p[x1l][x2-1][x3l]+nx3*nx4;
    p[x1l][x2][x3l][x4l]=p[x1l][x2-1][x3l][x4l]+nx3*nx4*nx5;
    p[x1l][x2][x3l][x4l][x5l]=p[x1l][x2-1][x3l][x4l][x5l]+nx3*nx4*nx5*nx6;
    for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2][x3l][x4l][x5]=p[x1l][x2][x3l][x4l][x5-1]+nx6;
    for(int x4=x4l+1;x4<=x4h;++x4){
      p[x1l][x2][x3l][x4]=p[x1l][x2][x3l][x4-1]+nx5;    
      p[x1l][x2][x3l][x4][x5l]=p[x1l][x2][x3l][x4-1][x5l]+nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2][x3l][x4][x5]=p[x1l][x2][x3l][x4][x5-1]+nx6;
    }
    for(int x3=x3l+1;x3<=x3h;++x3){
      p[x1l][x2][x3]=p[x1l][x2][x3-1]+nx4;
      p[x1l][x2][x3][x4l]=p[x1l][x2][x3-1][x4l]+nx4*nx5;
      p[x1l][x2][x3][x4l][x5l]=p[x1l][x2][x3-1][x4l][x5l]+nx4*nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2][x3][x4l][x5]=p[x1l][x2][x3][x4l][x5-1]+nx6;
      for(int x4=x4l+1;x4<=x4h;++x4){
        p[x1l][x2][x3][x4]=p[x1l][x2][x3][x4-1]+nx5;    
        p[x1l][x2][x3][x4][x5l]=p[x1l][x2][x3][x4-1][x5l]+nx5*nx6;
        for(int x5=x5l+1;x5<=x5h;++x5) p[x1l][x2][x3][x4][x5]=p[x1l][x2][x3][x4][x5-1]+nx6;
      }
    }
  }
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    p[x1][x2l][x3l]=p[x1-1][x2l][x3l]+nx2*nx3*nx4;
    p[x1][x2l][x3l][x4l]=p[x1-1][x2l][x3l][x4l]+nx2*nx3*nx4*nx5;
    p[x1][x2l][x3l][x4l][x5l]=p[x1-1][x2l][x3l][x4l][x5l]+nx2*nx3*nx4*nx5*nx6;
    for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2l][x3l][x4l][x5]=p[x1][x2l][x3l][x4l][x5-1]+nx6;
    for(int x4=x4l+1;x4<=x4h;++x4){
      p[x1][x2l][x3l][x4]=p[x1][x2l][x3l][x4-1]+nx5;    
      p[x1][x2l][x3l][x4][x5l]=p[x1][x2l][x3l][x4-1][x5l]+nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2l][x3l][x4][x5]=p[x1][x2l][x3l][x4][x5-1]+nx6;
    }
    for(int x3=x3l+1;x3<=x3h;++x3){
      p[x1][x2l][x3]=p[x1][x2l][x3-1]+nx4;
      p[x1][x2l][x3][x4l]=p[x1][x2l][x3-1][x4l]+nx4*nx5;
      p[x1][x2l][x3][x4l][x5l]=p[x1][x2l][x3-1][x4l][x5l]+nx4*nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2l][x3][x4l][x5]=p[x1][x2l][x3][x4l][x5-1]+nx6;
      for(int x4=x4l+1;x4<=x4h;++x4){
        p[x1][x2l][x3][x4]=p[x1][x2l][x3][x4-1]+nx5;    
        p[x1][x2l][x3][x4][x5l]=p[x1][x2l][x3][x4-1][x5l]+nx5*nx6;
        for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2l][x3][x4][x5]=p[x1][x2l][x3][x4][x5-1]+nx6;
      }
    }
    for(int x2=x2l+1;x2<=x2h;++x2){
      p[x1][x2]=p[x1][x2-1]+nx3;
      p[x1][x2][x3l]=p[x1][x2-1][x3l]+nx3*nx4;
      p[x1][x2][x3l][x4l]=p[x1][x2-1][x3l][x4l]+nx3*nx4*nx5;
      p[x1][x2][x3l][x4l][x5l]=p[x1][x2-1][x3l][x4l][x5l]+nx3*nx4*nx5*nx6;
      for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2][x3l][x4l][x5]=p[x1][x2][x3l][x4l][x5-1]+nx6;
      for(int x4=x4l+1;x4<=x4h;++x4){
        p[x1][x2][x3l][x4]=p[x1][x2][x3l][x4-1]+nx5;    
        p[x1][x2][x3l][x4][x5l]=p[x1][x2][x3l][x4-1][x5l]+nx5*nx6;
        for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2][x3l][x4][x5]=p[x1][x2][x3l][x4][x5-1]+nx6;
      }
      for(int x3=x3l+1;x3<=x3h;++x3){
        p[x1][x2][x3]=p[x1][x2][x3-1]+nx4;
        p[x1][x2][x3][x4l]=p[x1][x2][x3-1][x4l]+nx4*nx5;
        p[x1][x2][x3][x4l][x5l]=p[x1][x2][x3-1][x4l][x5l]+nx4*nx5*nx6;
        for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2][x3][x4l][x5]=p[x1][x2][x3][x4l][x5-1]+nx6;
        for(int x4=x4l+1;x4<=x4h;++x4){
          p[x1][x2][x3][x4]=p[x1][x2][x3][x4-1]+nx5;    
          p[x1][x2][x3][x4][x5l]=p[x1][x2][x3][x4-1][x5l]+nx5*nx6;
          for(int x5=x5l+1;x5<=x5h;++x5) p[x1][x2][x3][x4][x5]=p[x1][x2][x3][x4][x5-1]+nx6;
        }
      }
    }
  }
  return p;
}

void del_d6dim(double ******p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,int x6l,int x6h)
{
  delete[] (p[x1l][x2l][x3l][x4l][x5l]+x6l);
  delete[] (p[x1l][x2l][x3l][x4l]+x5l);
  delete[] (p[x1l][x2l][x3l]+x4l);
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

double *******d7dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,int x6l,int x6h,int x7l,int x7h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1,nx4=x4h-x4l+1,nx5=x5h-x5l+1,nx6=x6h-x6l+1,nx7=x7h-x7l+1;

  double *******p;
  p=new double****** [nx1] - x1l;
  p[x1l]=new double***** [nx1*nx2] - x2l;
  p[x1l][x2l]=new double**** [nx1*nx2*nx3] - x3l;
  p[x1l][x2l][x3l]=new double*** [nx1*nx2*nx3*nx4] - x4l;
  p[x1l][x2l][x3l][x4l]=new double** [nx1*nx2*nx3*nx4*nx5] - x5l;
  p[x1l][x2l][x3l][x4l][x5l]=new double* [nx1*nx2*nx3*nx4*nx5*nx6] - x6l;
  p[x1l][x2l][x3l][x4l][x5l][x6l]=new double [nx1*nx2*nx3*nx4*nx5*nx6*nx7] - x7l;

  for(int x6=x6l+1;x6<=x6h;++x6) p[x1l][x2l][x3l][x4l][x5l][x6]=p[x1l][x2l][x3l][x4l][x5l][x6-1]+nx7;
  for(int x5=x5l+1;x5<=x5h;++x5){
    p[x1l][x2l][x3l][x4l][x5]=p[x1l][x2l][x3l][x4l][x5-1]+nx6;    
    p[x1l][x2l][x3l][x4l][x5][x6l]=p[x1l][x2l][x3l][x4l][x5-1][x6l]+nx6*nx7;
    for(int x6=x6l+1;x6<=x6h;++x6) p[x1l][x2l][x3l][x4l][x5][x6]=p[x1l][x2l][x3l][x4l][x5][x6-1]+nx7;
  }
  for(int x4=x4l+1;x4<=x4h;++x4){
    p[x1l][x2l][x3l][x4]=p[x1l][x2l][x3l][x4-1]+nx5;    
    p[x1l][x2l][x3l][x4][x5l]=p[x1l][x2l][x3l][x4-1][x5l]+nx5*nx6;
    p[x1l][x2l][x3l][x4][x5l][x6l]=p[x1l][x2l][x3l][x4-1][x5l][x6l]+nx5*nx6*nx7;
    for(int x6=x6l+1;x6<=x6h;++x6) p[x1l][x2l][x3l][x4][x5l][x6]=p[x1l][x2l][x3l][x4][x5l][x6-1]+nx7;
    for(int x5=x5l+1;x5<=x5h;++x5){
    p[x1l][x2l][x3l][x4][x5]=p[x1l][x2l][x3l][x4][x5-1]+nx6;
    p[x1l][x2l][x3l][x4][x5][x6l]=p[x1l][x2l][x3l][x4][x5-1][x6l]+nx6*nx7;
    for(int x6=x6l+1;x6<=x6h;++x6) p[x1l][x2l][x3l][x4][x5][x6]=p[x1l][x2l][x3l][x4][x5][x6-1]+nx7;
   }
  }
  for(int x3=x3l+1;x3<=x3h;++x3){
    p[x1l][x2l][x3]=p[x1l][x2l][x3-1]+nx4;
    p[x1l][x2l][x3][x4l]=p[x1l][x2l][x3-1][x4l]+nx4*nx5;
    p[x1l][x2l][x3][x4l][x5l]=p[x1l][x2l][x3-1][x4l][x5l]+nx4*nx5*nx6;
    p[x1l][x2l][x3][x4l][x5l][x6l]=p[x1l][x2l][x3-1][x4l][x5l][x6l]+nx4*nx5*nx6*nx7;
    for(int x6=x6l+1;x6<=x6h;++x6) p[x1l][x2l][x3][x4l][x5l][x6]=p[x1l][x2l][x3][x4l][x5l][x6-1]+nx7;
    for(int x5=x5l+1;x5<=x5h;++x5){
      p[x1l][x2l][x3][x4l][x5]=p[x1l][x2l][x3][x4l][x5-1]+nx6;    
      p[x1l][x2l][x3][x4l][x5][x6l]=p[x1l][x2l][x3][x4l][x5-1][x6l]+nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1l][x2l][x3][x4l][x5][x6]=p[x1l][x2l][x3][x4l][x5][x6-1]+nx7;
    }
    for(int x4=x4l+1;x4<=x4h;++x4){
      p[x1l][x2l][x3][x4]=p[x1l][x2l][x3][x4-1]+nx5; 
      p[x1l][x2l][x3][x4][x5l]=p[x1l][x2l][x3][x4-1][x5l]+nx5*nx6;
      p[x1l][x2l][x3][x4][x5l][x6l]=p[x1l][x2l][x3][x4-1][x5l][x6l]+nx5*nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1l][x2l][x3][x4][x5l][x6]=p[x1l][x2l][x3][x4][x5l][x6-1]+nx7;
    for(int x5=x5l+1;x5<=x5h;++x5){
      p[x1l][x2l][x3][x4][x5]=p[x1l][x2l][x3][x4][x5-1]+nx6;
      p[x1l][x2l][x3][x4][x5][x6l]=p[x1l][x2l][x3][x4][x5-1][x6l]+nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1l][x2l][x3][x4][x5][x6]=p[x1l][x2l][x3][x4][x5][x6-1]+nx7;
     }
    }
  }

  for(int x2=x2l+1;x2<=x2h;++x2){
    p[x1l][x2]=p[x1l][x2-1]+nx3;
    p[x1l][x2][x3l]=p[x1l][x2-1][x3l]+nx3*nx4;
    p[x1l][x2][x3l][x4l]=p[x1l][x2-1][x3l][x4l]+nx3*nx4*nx5;
    p[x1l][x2][x3l][x4l][x5l]=p[x1l][x2-1][x3l][x4l][x5l]+nx3*nx4*nx5*nx6;
    p[x1l][x2][x3l][x4l][x5l][x6l]=p[x1l][x2-1][x3l][x4l][x5l][x6l]+nx3*nx4*nx5*nx6*nx7;
    for(int x6=x6l+1;x6<=x6h;++x6) p[x1l][x2][x3l][x4l][x5l][x6]=p[x1l][x2][x3l][x4l][x5l][x6-1]+nx7;
    for(int x5=x5l+1;x5<=x5h;++x5){
      p[x1l][x2][x3l][x4l][x5]=p[x1l][x2][x3l][x4l][x5-1]+nx6;
      p[x1l][x2][x3l][x4l][x5][x6l]=p[x1l][x2][x3l][x4l][x5-1][x6l]+nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1l][x2][x3l][x4l][x5][x6]=p[x1l][x2][x3l][x4l][x5][x6-1]+nx7;
     }
    for(int x4=x4l+1;x4<=x4h;++x4){
      p[x1l][x2][x3l][x4]=p[x1l][x2][x3l][x4-1]+nx5;    
      p[x1l][x2][x3l][x4][x5l]=p[x1l][x2][x3l][x4-1][x5l]+nx5*nx6;
      p[x1l][x2][x3l][x4][x5l][x6l]=p[x1l][x2][x3l][x4-1][x5l][x6l]+nx5*nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1l][x2][x3l][x4][x5l][x6]=p[x1l][x2][x3l][x4][x5l][x6-1]+nx7;
    }
    for(int x3=x3l+1;x3<=x3h;++x3){
      p[x1l][x2][x3]=p[x1l][x2][x3-1]+nx4;
      p[x1l][x2][x3][x4l]=p[x1l][x2][x3-1][x4l]+nx4*nx5;
      p[x1l][x2][x3][x4l][x5l]=p[x1l][x2][x3-1][x4l][x5l]+nx4*nx5*nx6;
      p[x1l][x2][x3][x4l][x5l][x6l]=p[x1l][x2][x3-1][x4l][x5l][x6l]+nx4*nx5*nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1l][x2][x3][x4l][x5l][x6]=p[x1l][x2][x3][x4l][x5l][x6-1]+nx7;
       }
      for(int x4=x4l+1;x4<=x4h;++x4){
        p[x1l][x2][x3l][x4]=p[x1l][x2][x3l][x4-1]+nx5;    
        p[x1l][x2][x3l][x4][x5l]=p[x1l][x2][x3l][x4-1][x5l]+nx5*nx6;
        p[x1l][x2][x3l][x4][x5l][x6l]=p[x1l][x2][x3l][x4-1][x5l][x6l]+nx5*nx6*nx7;
      for(int x5=x5l+1;x5<=x5h;++x5){
      p[x1l][x2][x3l][x4][x5]=p[x1l][x2][x3l][x4][x5-1]+nx6;
      p[x1l][x2][x3l][x4][x5][x6l]=p[x1l][x2][x3l][x4][x5-1][x6l]+nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1l][x2][x3l][x4][x5][x6]=p[x1l][x2][x3l][x4][x5][x6-1]+nx7;
       }
      }
    for(int x3=x3l+1;x3<=x3h;++x3){
     p[x1l][x2][x3]=p[x1l][x2][x3-1]+nx4;
     p[x1l][x2][x3][x4l]=p[x1l][x2][x3-1][x4l]+nx4*nx5;
     p[x1l][x2][x3][x4l][x5l]=p[x1l][x2][x3-1][x4l][x5l]+nx4*nx5*nx6;
     p[x1l][x2][x3][x4l][x5l][x6l]=p[x1l][x2][x3-1][x4l][x5l][x6l]+nx4*nx5*nx6*nx7;
    for(int x5=x5l+1;x5<=x5h;++x5){
      p[x1l][x2][x3][x4l][x5]=p[x1l][x2][x3][x4l][x5-1]+nx6;
      p[x1l][x2][x3][x4l][x5][x6l]=p[x1l][x2][x3][x4l][x5-1][x6l]+nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1l][x2][x3][x4l][x5][x6]=p[x1l][x2][x3][x4l][x5][x6-1]+nx7;
    }
    for(int x4=x4l+1;x4<=x4h;++x4){
     p[x1l][x2][x3][x4]=p[x1l][x2][x3][x4-1]+nx5;
     p[x1l][x2][x3][x4][x5l]=p[x1l][x2][x3][x4-1][x5l]+nx5*nx6;
     p[x1l][x2][x3][x4][x5l][x6l]=p[x1l][x2][x3][x4-1][x5l][x6l]+nx5*nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1l][x2][x3][x4][x5l][x6]=p[x1l][x2][x3][x4][x5l][x6-1]+nx7;
    }
    for(int x4=x4l+1;x4<=x4h;++x4){
     p[x1l][x2][x3][x4]=p[x1l][x2][x3][x4-1]+nx5;
     p[x1l][x2][x3][x4][x5l]=p[x1l][x2][x3][x4-1][x5l]+nx5*nx6;
     p[x1l][x2][x3][x4][x5l][x6l]=p[x1l][x2][x3][x4-1][x5l][x6l]+nx5*nx6*nx7;
      for(int x5=x5l+1;x5<=x5h;++x5){
      p[x1l][x2][x3][x4][x5]=p[x1l][x2][x3][x4][x5-1]+nx6;
      p[x1l][x2][x3][x4][x5][x6l]=p[x1l][x2][x3][x4][x5-1][x6l]+nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1l][x2][x3][x4][x5][x6]=p[x1l][x2][x3][x4][x5][x6-1]+nx7;
       }
      }
     }
    }

  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    p[x1][x2l][x3l]=p[x1-1][x2l][x3l]+nx2*nx3*nx4;
    p[x1][x2l][x3l][x4l]=p[x1-1][x2l][x3l][x4l]+nx2*nx3*nx4*nx5;
    p[x1][x2l][x3l][x4l][x5l]=p[x1-1][x2l][x3l][x4l][x5l]+nx2*nx3*nx4*nx5*nx6;
    p[x1][x2l][x3l][x4l][x5l][x6l]=p[x1-1][x2l][x3l][x4l][x5l][x6l]+nx2*nx3*nx4*nx5*nx6*nx7;
    for(int x6=x6l+1;x6<=x6h;++x6) p[x1][x2l][x3l][x4l][x5l][x6]=p[x1][x2l][x3l][x4l][x5l][x6-1]+nx7;
      for(int x5=x5l+1;x5<=x5h;++x5){
      p[x1][x2l][x3l][x4l][x5]=p[x1][x2l][x3l][x4l][x5-1]+nx6;
      p[x1][x2l][x3l][x4l][x5][x6l]=p[x1][x2l][x3l][x4l][x5-1][x6l]+nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1][x2l][x3l][x4l][x5][x6]=p[x1][x2l][x3l][x4l][x5][x6-1]+nx7;
       }
    for(int x4=x4l+1;x4<=x4h;++x4){
      p[x1][x2l][x3l][x4]=p[x1][x2l][x3l][x4-1]+nx5;
      p[x1][x2l][x3l][x4][x5l]=p[x1][x2l][x3l][x4-1][x5l]+nx5*nx6;
      p[x1][x2l][x3l][x4][x5l][x6l]=p[x1][x2l][x3l][x4-1][x5l][x6l]+nx5*nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1][x2l][x3l][x4][x5l][x6]=p[x1][x2l][x3l][x4][x5l][x6-1]+nx7;
    }
    for(int x3=x3l+1;x3<=x3h;++x3){
      p[x1][x2l][x3]=p[x1][x2l][x3-1]+nx4;
      p[x1][x2l][x3][x4l]=p[x1][x2l][x3-1][x4l]+nx4*nx5;
      p[x1][x2l][x3][x4l][x5l]=p[x1][x2l][x3-1][x4l][x5l]+nx4*nx5*nx6;
      p[x1][x2l][x3][x4l][x5l][x6l]=p[x1][x2l][x3-1][x4l][x5l][x6l]+nx4*nx5*nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1][x2l][x3][x4l][x5l][x6]=p[x1][x2l][x3][x4l][x5l][x6-1]+nx7;
       }
    for(int x2=x2l+1;x2<=x2h;++x2){
      p[x1][x2]=p[x1][x2-1]+nx3;
      p[x1][x2][x3l]=p[x1][x2-1][x3l]+nx3*nx4;
      p[x1][x2][x3l][x4l]=p[x1][x2-1][x3l][x4l]+nx3*nx4*nx5;
      p[x1][x2][x3l][x4l][x5l]=p[x1][x2-1][x3l][x4l][x5l]+nx3*nx4*nx5*nx6;
      p[x1][x2][x3l][x4l][x5l][x6l]=p[x1][x2-1][x3l][x4l][x5l][x6l]+nx3*nx4*nx5*nx6*nx7;
       for(int x6=x6l+1;x6<=x6h;++x6) p[x1][x2][x3l][x4l][x5l][x6]=p[x1][x2][x3l][x4l][x5l][x6-1]+nx7;
      }
      for(int x4=x4l+1;x4<=x4h;++x4){
        p[x1][x2l][x3l][x4]=p[x1][x2l][x3l][x4-1]+nx5;
        p[x1][x2l][x3l][x4][x5l]=p[x1][x2l][x3l][x4-1][x5l]+nx5*nx6;
        p[x1][x2l][x3l][x4][x5l][x6l]=p[x1][x2l][x3l][x4-1][x5l][x6l]+nx5*nx6*nx7;
      for(int x5=x5l+1;x5<=x5h;++x5){
      p[x1][x2l][x3l][x4][x5]=p[x1][x2l][x3l][x4][x5-1]+nx6;
      p[x1][x2l][x3l][x4][x5][x6l]=p[x1][x2l][x3l][x4][x5-1][x6l]+nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1][x2l][x3l][x4][x5][x6]=p[x1][x2l][x3l][x4][x5][x6-1]+nx7;
       }
      }
    for(int x3=x3l+1;x3<=x3h;++x3){
     p[x1][x2l][x3]=p[x1][x2l][x3-1]+nx4;
     p[x1][x2l][x3][x4l]=p[x1][x2l][x3-1][x4l]+nx4*nx5;
     p[x1][x2l][x3][x4l][x5l]=p[x1][x2l][x3-1][x4l][x5l]+nx4*nx5*nx6;
     p[x1][x2l][x3][x4l][x5l][x6l]=p[x1][x2l][x3-1][x4l][x5l][x6l]+nx4*nx5*nx6*nx7;
    for(int x5=x5l+1;x5<=x5h;++x5){
      p[x1][x2l][x3][x4l][x5]=p[x1][x2l][x3][x4l][x5-1]+nx6;
      p[x1][x2l][x3][x4l][x5][x6l]=p[x1][x2l][x3][x4l][x5-1][x6l]+nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1][x2l][x3][x4l][x5][x6]=p[x1][x2l][x3][x4l][x5][x6-1]+nx7;
    }
    for(int x4=x4l+1;x4<=x4h;++x4){
     p[x1][x2l][x3][x4]=p[x1][x2l][x3][x4-1]+nx5;
     p[x1][x2l][x3][x4][x5l]=p[x1][x2l][x3][x4-1][x5l]+nx5*nx6;
     p[x1][x2l][x3][x4][x5l][x6l]=p[x1][x2l][x3][x4-1][x5l][x6l]+nx5*nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1][x2l][x3][x4][x5l][x6]=p[x1][x2l][x3][x4][x5l][x6-1]+nx7;
     for(int x5=x5l+1;x5<=x5h;++x5){
      p[x1][x2l][x3][x4][x5]=p[x1][x2l][x3][x4][x5-1]+nx6;
      p[x1][x2l][x3][x4][x5][x6l]=p[x1][x2l][x3][x4][x5-1][x6l]+nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1][x2l][x3][x4][x5][x6]=p[x1][x2l][x3][x4][x5][x6-1]+nx7;
      }
     }
    }
    for(int x2=x2l+1;x2<=x2h;++x2){
     p[x1][x2][x3l]=p[x1][x2-1][x3l]+nx4;
     p[x1][x2][x3l][x4l]=p[x1][x2-1][x3l][x4l]+nx4*nx5;
     p[x1][x2][x3l][x4l][x5l]=p[x1][x2-1][x3l][x4l][x5l]+nx4*nx5*nx6;
     p[x1][x2][x3l][x4l][x5l][x6l]=p[x1][x2-1][x3l][x4l][x5l][x6l]+nx4*nx5*nx6*nx7;
    for(int x5=x5l+1;x5<=x5h;++x5){
      p[x1][x2][x3l][x4l][x5]=p[x1][x2][x3l][x4l][x5-1]+nx6;
      p[x1][x2][x3l][x4l][x5][x6l]=p[x1][x2][x3l][x4l][x5-1][x6l]+nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1][x2][x3l][x4l][x5][x6]=p[x1][x2][x3l][x4l][x5][x6-1]+nx7;
    }
    for(int x4=x4l+1;x4<=x4h;++x4){
     p[x1][x2][x3l][x4]=p[x1][x2][x3l][x4-1]+nx5;
     p[x1][x2][x3l][x4][x5l]=p[x1][x2][x3l][x4-1][x5l]+nx5*nx6;
     p[x1][x2][x3l][x4][x5l][x6l]=p[x1][x2][x3l][x4-1][x5l][x6l]+nx5*nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1][x2][x3l][x4][x5l][x6]=p[x1][x2][x3l][x4][x5l][x6-1]+nx7;
     for(int x5=x5l+1;x5<=x5h;++x5){
      p[x1][x2][x3l][x4][x5]=p[x1][x2][x3l][x4][x5-1]+nx6;
      p[x1][x2][x3l][x4][x5][x6l]=p[x1][x2][x3l][x4][x5-1][x6l]+nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1][x2][x3l][x4][x5][x6]=p[x1][x2][x3l][x4][x5][x6-1]+nx7;
     }
    }
    for(int x3=x3l+1;x3<=x3h;++x3){
      p[x1][x2][x3]=p[x1][x2][x3-1]+nx4;
      p[x1][x2][x3][x4l]=p[x1][x2][x3-1][x4l]+nx4*nx5;
      p[x1][x2][x3][x4l][x5l]=p[x1][x2][x3-1][x4l][x5l]+nx4*nx5*nx6;
      p[x1][x2][x3][x4l][x5l][x6l]=p[x1][x2][x3-1][x4l][x5l][x6l]+nx4*nx5*nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1][x2][x3][x4l][x5l][x6]=p[x1][x2][x3][x4l][x5l][x6-1]+nx7;
     for(int x5=x5l+1;x5<=x5h;++x5){
      p[x1][x2][x3][x4l][x5]=p[x1][x2][x3][x4l][x5-1]+nx6;
      p[x1][x2][x3][x4l][x5][x6l]=p[x1][x2][x3][x4l][x5-1][x6l]+nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1][x2][x3][x4l][x5][x6]=p[x1][x2][x3][x4l][x5][x6-1]+nx7;
      }
    for(int x4=x4l+1;x4<=x4h;++x4){
     p[x1][x2][x3][x4]=p[x1][x2][x3][x4-1]+nx5;
     p[x1][x2][x3][x4][x5l]=p[x1][x2][x3][x4-1][x5l]+nx5*nx6;
     p[x1][x2][x3][x4][x5l][x6l]=p[x1][x2][x3][x4-1][x5l][x6l]+nx5*nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1][x2][x3][x4][x5l][x6]=p[x1][x2][x3][x4][x5l][x6-1]+nx7;
     for(int x5=x5l+1;x5<=x5h;++x5){
      p[x1][x2][x3][x4][x5]=p[x1][x2][x3][x4][x5-1]+nx6;
      p[x1][x2][x3][x4][x5][x6l]=p[x1][x2][x3][x4][x5-1][x6l]+nx6*nx7;
      for(int x6=x6l+1;x6<=x6h;++x6) p[x1][x2][x3][x4][x5][x6]=p[x1][x2][x3][x4][x5][x6-1]+nx7;
      }
     }
    }
   }
  }
  return p;
}

void del_d7dim(double *******p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l,int x5h,int x6l,int x6h,
int x7l,int x7h)
{
  delete[] (p[x1l][x2l][x3l][x4l][x5l][x6l]+x7l);
  delete[] (p[x1l][x2l][x3l][x4l][x5l]+x6l);
  delete[] (p[x1l][x2l][x3l][x4l]+x5l);
  delete[] (p[x1l][x2l][x3l]+x4l);
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

int *i1dim(int x1l,int x1h)
{
  return new int [x1h-x1l+1] - x1l;
}

void del_i1dim(int *p,int x1l,int x1h)
{
  delete[] (p+x1l);
}

int **i2dim(int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  int **p;
  p=new int* [nx1] - x1l;
  p[x1l]=new int [nx1*nx2] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

void del_i2dim(int **p,int x1l, int x1h, int x2l, int x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

int ***i3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  int ***p;
  p=new int** [nx1] - x1l;
  p[x1l]=new int* [nx1*nx2] - x2l;
  p[x1l][x2l]=new int [nx1*nx2*nx3] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

void del_i3dim(int ***p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

int ***** i5dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l, int x5h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1,nx4=x4h-x4l+1,nx5=x5h-x5l+1;
  int *****p;
  p=new int**** [nx1] - x1l;
  p[x1l]=new int*** [nx1*nx2] - x2l;
  p[x1l][x2l]=new int** [nx1*nx2*nx3] - x3l;
  p[x1l][x2l][x3l]=new int* [nx1*nx2*nx3*nx4] - x4l;
  p[x1l][x2l][x3l][x4l]=new int [nx1*nx2*nx3*nx4*nx5] - x5l;
  
  for(int x4=x4l+1;x4<=x4h;++x4) p[x1l][x2l][x3l][x4]=p[x1l][x2l][x3l][x4-1]+nx5;
  for(int x3=x3l+1;x3<=x3h;++x3){
    p[x1l][x2l][x3]=p[x1l][x2l][x3-1]+nx4;    
    p[x1l][x2l][x3][x4l]=p[x1l][x2l][x3-1][x4l]+nx4*nx5;
    for(int x4=x4l+1;x4<=x4h;++x4) p[x1l][x2l][x3][x4]=p[x1l][x2l][x3][x4-1]+nx5;
  }
  for(int x2=x2l+1;x2<=x2h;++x2){
    p[x1l][x2]=p[x1l][x2-1]+nx3;
    p[x1l][x2][x3l]=p[x1l][x2-1][x3l]+nx3*nx4;
    p[x1l][x2][x3l][x4l]=p[x1l][x2-1][x3l][x4l]+nx3*nx4*nx5;
    for(int x4=x4l+1;x4<=x4h;++x4) p[x1l][x2][x3l][x4]=p[x1l][x2][x3l][x4-1]+nx5;
    for(int x3=x3l+1;x3<=x3h;++x3){
      p[x1l][x2][x3]=p[x1l][x2][x3-1]+nx4;
      p[x1l][x2][x3][x4l]=p[x1l][x2][x3-1][x4l]+nx4*nx5;
      for(int x4=x4l+1;x4<=x4h;++x4) p[x1l][x2][x3][x4]=p[x1l][x2][x3][x4-1]+nx5;
    }
  }	
  for(int x1=x1l+1;x1<=x1h;++x1) {
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    p[x1][x2l][x3l]=p[x1-1][x2l][x3l]+nx2*nx3*nx4;
    p[x1][x2l][x3l][x4l]=p[x1-1][x2l][x3l][x4l]+nx2*nx3*nx4*nx5;
    for(int x4=x4l+1;x4<=x4h;++x4) p[x1][x2l][x3l][x4]=p[x1][x2l][x3l][x4-1]+nx5;
    for(int x3=x3l+1;x3<=x3h;++x3){
      p[x1][x2l][x3]=p[x1][x2l][x3-1]+nx4;
      p[x1][x2l][x3][x4l]=p[x1][x2l][x3-1][x4l]+nx4*nx5;
      for(int x4=x4l+1;x4<=x4h;++x4) p[x1][x2l][x3][x4]=p[x1][x2l][x3][x4-1]+nx5;
    }
    for(int x2=x2l+1;x2<=x2h;++x2){
      p[x1][x2]=p[x1][x2-1]+nx3;
      p[x1][x2][x3l]=p[x1][x2-1][x3l]+nx3*nx4;
      p[x1][x2][x3l][x4l]=p[x1][x2-1][x3l][x4l]+nx3*nx4*nx5;
      for(int x4=x4l+1;x4<=x4h;++x4) p[x1][x2][x3l][x4]=p[x1][x2][x3l][x4-1]+nx5;
      for(int x3=x3l+1;x3<=x3h;++x3){
        p[x1][x2][x3]=p[x1][x2][x3-1]+nx4;
        p[x1][x2][x3][x4l]=p[x1][x2][x3-1][x4l]+nx4*nx5;
        for(int x4=x4l+1;x4<=x4h;++x4) p[x1][x2][x3][x4]=p[x1][x2][x3][x4-1]+nx5;
      }
    }
  }
  return p;
}

void del_i5dim(int *****p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h,int x4l,int x4h,int x5l, int x5h)
{
  delete[] (p[x1l][x2l][x3l][x4l]+x5l);
  delete[] (p[x1l][x2l][x3l]+x4l);
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

int **i2dim(int x1l,int x1h,int x2l,int *x2h)
{
  int nx1=x1h-x1l+1,nx12=0;
  int **p;
  p=new int* [nx1] - x1l;
  for(int x1=x1l;x1<=x1h;++x1) nx12+=(x2h[x1]-x2l+1);
  p[x1l]=new int [nx12] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+(x2h[x1-1]-x2l+1);
  return p;
}

void del_i2dim(int **p,int x1l,int x1h,int x2l,int *x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

int ***i3dim(int x1l,int x1h,int x2l,int *x2h,int x3l,int **x3h)
{
  int nx1=x1h-x1l+1;
  int *nx2=new int [nx1]-x1l,nx12=0;
  int **nx3=i2dim(x1l,x1h,x2l,x2h),nx123=0;
  for(int x1=x1l;x1<=x1h;++x1){
    nx12+=(nx2[x1]=x2h[x1]-x2l+1);
    for(int x2=x2l;x2<=x2h[x1];++x2) nx123+=(nx3[x1][x2]=x3h[x1][x2]-x3l+1);
  }
//
  int ***p;
  p=new int** [nx1]-x1l;
  p[x1l]=new int* [nx12]-x2l;
  p[x1l][x2l]=new int [nx123]-x3l;
  for(int x2=x2l+1;x2<=x2h[x1l];++x2) p[x1l][x2]=p[x1l][x2-1]+nx3[x1l][x2-1];
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2[x1-1];
    p[x1][x2l]=p[x1-1][x2h[x1-1]]+nx3[x1-1][x2h[x1-1]];
    for(int x2=x2l+1;x2<=x2h[x1];++x2) p[x1][x2]=p[x1][x2-1]+nx3[x1][x2-1];
  }
//
  del_i2dim(nx3,x1l,x1h,x2l,x2h);
  delete[] (nx2+x1l);
//
  return p;
}

void del_i3dim(int ***p,int x1l,int x1h,int x2l,int *x2h,int x3l,int **x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

int ***i3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int *x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nxt=0;
  int ***p;
  p=new int** [nx1] - x1l;
  p[x1l]=new int* [nx1*nx2] - x2l;
  int *nx3=new int [nx1] -x1l;
  for(int x1=x1l;x1<=x1h;++x1) nxt+=(nx3[x1]=x3h[x1]-x3l+1);
  p[x1l][x2l]=new int [nx2*nxt] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3[x1l];
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3[x1-1];
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3[x1];
  }
  delete[] (nx3+x1l);
  return p;
}

void del_i3dim(int ***p,int x1l,int x1h,int x2l,int x2h,int x3l,int *x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

void ***v2dim(int x1l,int x1h,int x2l,int x2h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1;
  void ***p;
  p=new void** [nx1] - x1l;
  p[x1l]=new void* [nx1*nx2] - x2l;
  for(int x1=x1l+1;x1<=x1h;++x1) p[x1]=p[x1-1]+nx2;
  return p;
}

void del_v2dim(void ***p,int x1l, int x1h, int x2l, int x2h)
{
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

void ****v3dim(int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  int nx1=x1h-x1l+1,nx2=x2h-x2l+1,nx3=x3h-x3l+1;
  void ****p;
  p=new void*** [nx1] - x1l;
  p[x1l]=new void** [nx1*nx2] - x2l;
  p[x1l][x2l]=new void* [nx1*nx2*nx3] - x3l;
  for(int x2=x2l+1;x2<=x2h;++x2) p[x1l][x2]=p[x1l][x2-1]+nx3;
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2;
    p[x1][x2l]=p[x1-1][x2l]+nx2*nx3;
    for(int x2=x2l+1;x2<=x2h;++x2) p[x1][x2]=p[x1][x2-1]+nx3;
  }
  return p;
}

void del_v3dim(void ****p,int x1l,int x1h,int x2l,int x2h,int x3l,int x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}

void ****v3dim(int x1l,int x1h,int x2l,int *x2h,int x3l,int **x3h)
{
  int nx1=x1h-x1l+1;
  int *nx2=new int [nx1]-x1l,nx12=0;
  int **nx3=i2dim(x1l,x1h,x2l,x2h),nx123=0;
  for(int x1=x1l;x1<=x1h;++x1){
    nx12+=(nx2[x1]=x2h[x1]-x2l+1);
    for(int x2=x2l;x2<=x2h[x1];++x2) nx123+=(nx3[x1][x2]=x3h[x1][x2]-x3l+1);
  }
//
  void ****p;
  p=new void*** [nx1]-x1l;
  p[x1l]=new void** [nx12]-x2l;
  p[x1l][x2l]=new void* [nx123]-x3l;
  for(int x2=x2l+1;x2<=x2h[x1l];++x2) p[x1l][x2]=p[x1l][x2-1]+nx3[x1l][x2-1];
  for(int x1=x1l+1;x1<=x1h;++x1){
    p[x1]=p[x1-1]+nx2[x1-1];
    p[x1][x2l]=p[x1-1][x2h[x1-1]]+nx3[x1-1][x2h[x1-1]];
    for(int x2=x2l+1;x2<=x2h[x1];++x2) p[x1][x2]=p[x1][x2-1]+nx3[x1][x2-1];
  }
//
  del_i2dim(nx3,x1l,x1h,x2l,x2h);
  delete[] (nx2+x1l);
//
  return p;
}

void del_v3dim(void ****p,int x1l,int x1h,int x2l,int *x2h,int x3l,int **x3h)
{
  delete[] (p[x1l][x2l]+x3l);
  delete[] (p[x1l]+x2l);
  delete[] (p+x1l);
}
