#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include "dfparser.h"

using namespace std;

void align(char* str) {
  int i;
  while( isspace(str[0]) ) {
    for(i=0;i<(int)strlen(str);i++)
      str[i] = str[i+1];
  }
}

void rmeqsign(char* str) {
  int i;
  for(i=0;i<(int)strlen(str);i++)
    if( str[i] == '=' )
      str[i] = ' ';
}

void rmbrackets(char* str) {
  int i;
  for(i=0;i<(int)strlen(str);i++)
    if( str[i] == '[' || str[i] == ']' )
      str[i] = ' ';
}

void cutword(char* str, char* word) {
  int n;
  char* blank;
  blank = strchr(str,' ');
  if( blank==NULL )
    n = strlen(str);
  else
    n = (int)(blank-str);
  word[n] = '\0';
  while( --n>=0 ) {
    word[n] = str[n];
    str[n] = ' ';
  }
  align(str);
}

int gettc(char* vartype, const char* type) {
  int n;
  char tmpstr[16];
  strcpy(tmpstr,type);
  align(tmpstr);
  rmbrackets(tmpstr);
  cutword(tmpstr,vartype);
  n = atoi(tmpstr);
  return ( n>1 ) ? n : 1;
}

int getvar_s(void* value, const char* name,
	     const char* type, const char* filename) {
  FILE* fptr;
  int error;

  if( (fptr=fopen(filename,"r"))!=NULL ) {
    error = fgetvar(value,name,type,fptr);
    fclose(fptr);
    if(error) {
      cout << name << " not set in parameters.dat ... abort" << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    return error;
  }
  else
    return -1;
}

int getvar(void* value, const char* name,
	   const char* type, const char* filename) {
  FILE* fptr;
  int error;

  if( (fptr=fopen(filename,"r"))!=NULL ) {
    error = fgetvar(value,name,type,fptr);
    fclose(fptr);
    return error;
  }
  else
    return -1;
}

int fgetvar(void* value, const char* name,
	    const char* type, FILE* fptr) {
  int error = 1;
  int i,count;

  char varname[64];
  char vartype[16];
  char buf[256];

  rewind(fptr);

  while( fgets(buf,256,fptr)!=NULL ) {
    align(buf);
    if( strlen(buf)!=0 && buf[0]!='#' && buf[0]!='%' && buf[0]!='!' ) {
      rmeqsign(buf);
      cutword(buf,varname);
      if( strcmp(varname,name)==0 ) {
	error = 0;
	count = gettc(vartype,type);
	if( strcmp(vartype,"char")==0 )
	  if( count==1 ) {
	    if( sscanf(buf,"%c",((char*)value))==EOF )
	      error = 1;
	  }
	  else {
	    if( sscanf(buf,"%s",((char*)value))==EOF )
	      error = 1;
	  }
	else if( strcmp(vartype,"char*")==0 ) {
	  if( sscanf(buf,"%s",((char*)value))==EOF )
	    error = 1;
	}
	else if( strcmp(vartype,"int")==0 ) {
	  for(i=0;i<count;i++) {
	    if( sscanf(buf,"%d",((int*)value+(size_t)i))==EOF )
	      error = 1;
	    cutword(buf,vartype);
	  }
	}
	else if( strcmp(vartype,"uint")==0 ) {
	  for(i=0;i<count;i++) {
	    if( sscanf(buf,"%u",((unsigned int*)value+(size_t)i))==EOF )
	      error = 1;
	    cutword(buf,vartype);
	  }
	}
	else if( strcmp(vartype,"short")==0 ) {
	  for(i=0;i<count;i++) {
	    if( sscanf(buf,"%hd",((short*)value+(size_t)i))==EOF )
	      error = 1;
	    cutword(buf,vartype);
	  }
	}
	else if( strcmp(vartype,"ushort")==0 ) {
	  for(i=0;i<count;i++) {
	    if( sscanf(buf,"%hu",((unsigned short*)value+(size_t)i))==EOF )
	      error = 1;
	    cutword(buf,vartype);
	  }
	}
	else if( strcmp(vartype,"long")==0 ) {
	  for(i=0;i<count;i++) {
	    if( sscanf(buf,"%ld",((long*)value+(size_t)i))==EOF )
	      error = 1;
	    cutword(buf,vartype);
	  }
	}
	else if( strcmp(vartype,"ulong")==0 ) {
	  for(i=0;i<count;i++) {
	    if( sscanf(buf,"%lu",((unsigned long*)value+(size_t)i))==EOF )
	      error = 1;
	    cutword(buf,vartype);
	  }
	}
	else if( strcmp(vartype,"float")==0 ) {
	  for(i=0;i<count;i++) {
	    if( sscanf(buf,"%g",((float*)value+(size_t)i))==EOF )
	      error = 1;
	    cutword(buf,vartype);
	  }
	}
	else if( strcmp(vartype,"double")==0 ) {
	  for(i=0;i<count;i++) {
	    if( sscanf(buf,"%lg",((double*)value+(size_t)i))==EOF )
	      error = 1;
	    cutword(buf,vartype);
	  }
	}
      }
    }
  }
  return error;
}

int getvarlist(dfpvar* list, const int N, const char* filename) {
  int i,error = 0;
  FILE *fptr;

  if( (fptr=fopen(filename,"r"))!=NULL ) {
    for(i=0;i<N;i++)
      error += fgetvar(list[i].value,
		       list[i].name,
		       list[i].type,
		       fptr);
    fclose(fptr);
    return error;
  }
  else
    return -1;
}

void setdfpvar(dfpvar* var, void* value, const char* name,
	       const char* type) {
  var->value = value;
  strcpy(var->name,name);
  strcpy(var->type,type);
}
