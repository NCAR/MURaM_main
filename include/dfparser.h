#ifndef __DFPARSER_INCLUDED__
#define __DFPARSER_INCLUDED__

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  char name[64];
  char type[16];
  void* value;
} dfpvar;

extern int getvar(void*, const char*, const char*, const char*);
extern int getvar_s(void*, const char*, const char*, const char*);
extern int fgetvar(void*, const char*, const char*, FILE*);
extern int getvarlist(dfpvar*, const int, const char*);
extern void setdfpvar(dfpvar*, void*, const char*, const char*);

#ifdef __cplusplus
}
#endif

#endif
