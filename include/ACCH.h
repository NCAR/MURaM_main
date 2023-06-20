#ifndef __ACCH_HEADER__
#define __ACCH_HEADER__

// NVIDIA NV TOOLS Profiler Push/Pop
#ifdef USE_NVTX
#include "nvToolsExt.h"
const uint32_t ProfilerColors[] = {
  0xff00ff00, 0xff0000ff, 0xffffff00, 0xffff00ff, 0xff00ffff, 0xffff0000,
  0xffffffff
};
#define NVPROF_PUSH_RANGE(name,cid) { \
    int color_id = cid%7; \
    nvtxEventAttributes_t eventAttrib = {0}; \
    eventAttrib.version = NVTX_VERSION; \
    eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE; \
    eventAttrib.colorType = NVTX_COLOR_ARGB; \
    eventAttrib.color = ProfilerColors[color_id]; \
    eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII; \
    eventAttrib.message.ascii = name; \
    nvtxRangePushEx(&eventAttrib); \
}
#define NVPROF_POP_RANGE nvtxRangePop();
#else
#define NVPROF_PUSH_RANGE(name,cid)
#define NVPROF_POP_RANGE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <iostream>

namespace ACCH
{

#ifdef CACHE_SIZE
const int shared_size = CACHE_SIZE;
#else
extern int shared_size;
#endif

void * Malloc(std::size_t numbytes);
void Free(void * ptr, std::size_t numbytes);

void * GPUMalloc(std::size_t numbytes);
void GPUFree(void * ptr);
void GPUFree(void * ptr, std::size_t numbytes);

void UpdateCPU(void * ptr, std::size_t numbytes);
void UpdateGPU(void * ptr, std::size_t numbytes);
void UpdateCPU(void * ptr, std::size_t start, std::size_t numbytes);
void UpdateGPU(void * ptr, std::size_t start, std::size_t numbytes);

void MemcpyCPUtoGPU(void * h_ptr, void * d_ptr, std::size_t numbytes);
void MemcpyGPUtoCPU(void * h_ptr, void * d_ptr, std::size_t numbytes);

void MemcpyCPU(void * src, void * dest, std::size_t bytes);
void MemcpyGPU(void * src, void * dest, std::size_t bytes);

void Create(void * ptr, std::size_t bytes);
void Copyin(void * ptr, std::size_t bytes);
void Delete(void * ptr, std::size_t bytes);
void Copyout(void * ptr, std::size_t bytes);

void Compare(
  void * data, const char * datatype, std::size_t numelements,
  const char * dataname, const char * filename, 
  const char * funcname, unsigned int linenum
);
void GPUCompare(
  void * data, const char * datatype, std::size_t numelements,
  std::size_t numbytes,
  const char * dataname, const char * filename, 
  const char * funcname, unsigned int linenum
);

bool Present(void * ptr, std::size_t bytes);

int GetNumGPU();
void SetGPU(int num);
int GetGPU();

void * GetDevicePtr(void * h_ptr);

template <typename T>
T* Malloc1D(
  std::size_t D
)
{
  return Malloc(D*sizeof(T));
}

template <typename T>
void Free1D(
  T * ptr, std::size_t D
)
{
  Free(ptr, D*sizeof(T));
}

template <typename T>
T* GPUMalloc1D(
  std::size_t D
)
{
  return GPUMalloc(D*sizeof(T));
}

template <typename T>
void GPUFree1D(
  T * ptr, std::size_t D
)
{
  GPUFree(ptr, D*sizeof(T));
}

template <typename T>
void UpdateCPU1D(
  T * ptr, std::size_t size
)
{
  UpdateCPU(ptr, size*sizeof(T));
}

template <typename T>
void UpdateGPU1D(
  T * ptr, std::size_t size
)
{
  UpdateGPU(ptr, size*sizeof(T));
}

template <typename T>
void UpdateCPU1D(
  T * ptr, std::size_t ind1, std::size_t ind2
)
{
  UpdateCPU(ptr, ind1*sizeof(T), ind2*sizeof(T));
}

template <typename T>
void UpdateGPU1D(
  T * ptr, std::size_t ind1, std::size_t ind2
)
{
  UpdateGPU(ptr, ind1*sizeof(T), ind2*sizeof(T));
}

template <typename T>
void MemcpyCPUtoGPU1D(
  T * h_ptr, T * d_ptr, std::size_t numelements
)
{
  MemcpyCPUtoGPU(h_ptr, d_ptr, numelements*sizeof(T));
}

template <typename T>
void MemcpyGPUtoCPU1D(
  T * h_ptr, T * d_ptr, std::size_t numelements
)
{
  MemcpyGPUtoCPU(h_ptr, d_ptr, numelements*sizeof(T));
}

template <typename T>
T** Malloc2D(
  std::size_t D1, std::size_t D2
)
{
  T* ptr =  (T*) Malloc(D1*D2*sizeof(T));
  T** ptr2D = (T**) Malloc(D1*sizeof(T*));
  for(int i = 0; i < D1; i++) {
    ptr2D[i] = ptr + i*D2;
  }
#ifdef _OPENACC
  T* d_ptr = (T*) GetDevicePtr(ptr);
  T** d_ptr2D = (T**) GetDevicePtr(ptr2D);
  T** tmp = (T**) malloc(D1*sizeof(T*));
  for(int i = 0; i < D1; i++) {
    tmp[i] = d_ptr + i*D2;
  }
  MemcpyCPUtoGPU(tmp, d_ptr2D, D1*sizeof(T*));
  free(tmp);
#endif
  return ptr2D;
}

template <typename T>
void Free2D(
  T** ptr, std::size_t D1, std::size_t D2
)
{
  Free(ptr[0], D1*D2*sizeof(T));
  Free(ptr, D1*sizeof(T*));
}

template<typename T>
void UpdateCPU2D(
  T** ptr, std::size_t D1, std::size_t D2
)
{
  UpdateCPU(ptr[0], D1*D2*sizeof(T));
}

template <typename T>
void UpdateGPU2D(
  T** ptr, std::size_t D1, std::size_t D2
)
{
  UpdateGPU(ptr[0], D1*D2*sizeof(T));
}

template <typename T>
T*** Malloc3D(
  std::size_t D1, std::size_t D2, std::size_t D3
)
{
  T * ptr = (T*) Malloc(D1*D2*D3*sizeof(T));
  T** ptr2D = (T**) Malloc(D1*D2*sizeof(T*));
  T*** ptr3D = (T***) Malloc(D1*sizeof(T**));
  for(int i = 0; i < D1; i++)
    for(int j = 0; j < D2; j++)
      ptr2D[i*D2+j] = ptr + (i*D2+j)*D3;
  for(int i = 0; i < D1; i++)
    ptr3D[i] = ptr2D + i*D2;
#ifdef _OPENACC
  T * d_ptr = (T*) GetDevicePtr(ptr);
  T** d_ptr2D = (T**) GetDevicePtr(ptr2D);
  T*** d_ptr3D = (T***) GetDevicePtr(ptr3D);
  T** tmp = (T**) malloc(D1*D2*sizeof(T*));
  for(int i = 0; i < D1; i++)
    for(int j = 0; j < D2; j++)
      tmp[i*D2+j] = d_ptr + (i*D2+j)*D3;
  MemcpyCPUtoGPU(tmp, d_ptr2D, D1*D2*sizeof(T*));
  free(tmp);
  T*** tmp2 = (T***) malloc(D1*sizeof(T**));
  for(int i = 0; i < D1; i++)
    tmp2[i] = d_ptr2D + i*D2;
  MemcpyCPUtoGPU(tmp2, d_ptr3D, D1*sizeof(T**));
  free(tmp2);
#endif
  return ptr3D;
}

template <typename T>
void Free3D(
  T*** ptr, std::size_t D1, std::size_t D2, std::size_t D3
)
{
  Free(ptr[0][0], D1*D2*D3*sizeof(T));
  Free(ptr[0], D1*D2*sizeof(T*));
  Free(ptr, D1*sizeof(T**));
}

template <typename T>
void UpdateCPU3D(
  T*** ptr, std::size_t D1, std::size_t D2, std::size_t D3
)
{
  UpdateCPU(ptr[0][0], D1*D2*D3*sizeof(T));
}

template <typename T>
void UpdateGPU3D(
  T*** ptr, std::size_t D1, std::size_t D2, std::size_t D3
)
{
  UpdateGPU(ptr[0][0], D1*D2*D3*sizeof(T));
}

template <typename T>
T**** Malloc4D(
  std::size_t D1, std::size_t D2,
  std::size_t D3, std::size_t D4
)
{
  T * ptr = (T*) Malloc(D1*D2*D3*D4*sizeof(T));
  T** ptr2D = (T**) Malloc(D1*D2*D3*sizeof(T*));
  T*** ptr3D = (T***) Malloc(D1*D2*sizeof(T**));
  T**** ptr4D = (T****) Malloc(D1*sizeof(T***));
  for(int i = 0; i < D1; i++)
    for(int j = 0; j < D2; j++)
      for(int k = 0; k < D3; k++)
        ptr2D[i*D2*D3+j*D3+k] = ptr + (i*D2*D3+j*D3+k)*D4;
  for(int i = 0; i < D1; i++)
    for(int j = 0; j < D2; j++)
      ptr3D[i*D2+j] = ptr2D + (i*D2+j)*D3;
  for(int i = 0; i < D1; i++)
    ptr4D[i] = ptr3D + i*D2;
#ifdef _OPENACC
  T * d_ptr = (T*) GetDevicePtr(ptr);
  T** d_ptr2D = (T**) GetDevicePtr(ptr2D);
  T*** d_ptr3D = (T***) GetDevicePtr(ptr3D);
  T**** d_ptr4D = (T****) GetDevicePtr(ptr4D);
  T** tmp2D = (T**) malloc(D1*D2*D3*sizeof(T*));
  T*** tmp3D = (T***) malloc(D1*D2*sizeof(T**));
  T**** tmp4D = (T****) malloc(D1*sizeof(T***));
  for(int i = 0; i < D1; i++)
    for(int j = 0; j < D2; j++)
      for(int k = 0; k < D3; k++)
        tmp2D[i*D2*D3+j*D3+k] = d_ptr + (i*D2*D3+j*D3+k)*D4;
  MemcpyCPUtoGPU(tmp2D, d_ptr2D, D1*D2*D3*sizeof(T*));
  for(int i = 0; i < D1; i++)
    for(int j = 0; j < D2; j++)
      tmp3D[i*D2+j] = d_ptr2D + (i*D2+j)*D3;
  MemcpyCPUtoGPU(tmp3D, d_ptr3D, D1*D2*sizeof(T**));
  for(int i = 0; i < D1; i++)
    tmp4D[i] = d_ptr3D + i*D2;
  MemcpyCPUtoGPU(tmp4D, d_ptr4D, D1*sizeof(T***));
  free(tmp2D);
  free(tmp3D);
  free(tmp4D);
#endif
  return ptr4D;
}

template <typename T>
void Free4D(
  T**** ptr, std::size_t D1, std::size_t D2,
  std::size_t D3, std::size_t D4
)
{
  Free(ptr[0][0][0], D1*D2*D3*D4*sizeof(T));
  Free(ptr[0][0], D1*D2*D3*sizeof(T*));
  Free(ptr[0], D1*D2*sizeof(T**));
  Free(ptr, D1*sizeof(T***));
}

template <typename T>
void UpdateCPU4D(
  T**** ptr, std::size_t D1, std::size_t D2, std::size_t D3,
  std::size_t D4
)
{
  UpdateCPU(ptr[0][0][0], D1*D2*D3*D4*sizeof(T));
}

template <typename T>
void UpdateGPU4D(
  T**** ptr, std::size_t D1, std::size_t D2, std::size_t D3,
  std::size_t D4
)
{
  UpdateGPU(ptr[0][0][0], D1*D2*D3*D4*sizeof(T));
}


template <typename T>
T ****** Malloc6D(
  std::size_t D1, std::size_t D2, std::size_t D3,
  std::size_t D4, std::size_t D5, std::size_t D6
)
{
  T * ptr = (T*) Malloc(D1*D2*D3*D4*D5*D6*sizeof(T));
  T ** ptr2D = (T**) Malloc(D1*D2*D3*D4*D5*sizeof(T*));
  T *** ptr3D = (T***) Malloc(D1*D2*D3*D4*sizeof(T**));
  T **** ptr4D = (T****) Malloc(D1*D2*D3*sizeof(T***));
  T ***** ptr5D = (T*****) Malloc(D1*D2*sizeof(T****));
  T ****** ptr6D = (T******) Malloc(D1*sizeof(T*****));

  for(int d1 = 0; d1 < D1; d1++)
    for(int d2 = 0; d2 < D2; d2++)
      for(int d3 = 0; d3 < D3; d3++)
        for(int d4 = 0; d4 < D4; d4++)
          for(int d5 = 0; d5 < D5; d5++) {
            int ind = d1*D2*D3*D4*D5 +
                      d2*   D3*D4*D5 +
                      d3*      D4*D5 +
                      d4*         D5 +
                      d5;
            ptr2D[ind] = ptr + ind*D6;
          }
  for(int d1 = 0; d1 < D1; d1++)
    for(int d2 = 0; d2 < D2; d2++)
      for(int d3 = 0; d3 < D3; d3++)
        for(int d4 = 0; d4 < D4; d4++) {
          int ind = d1*D2*D3*D4 +
                    d2*   D3*D4 +
                    d3*      D4 +
                    d4;
          ptr3D[ind] = ptr2D + ind*D5;
        }
  for(int d1 = 0; d1 < D1; d1++)
    for(int d2 = 0; d2 < D2; d2++)
      for(int d3 = 0; d3 < D3; d3++) {
        int ind = d1*D2*D3 +
                  d2*   D3 +
                  d3;
        ptr4D[ind] = ptr3D + ind*D4;
      }
  for(int d1 = 0; d1 < D1; d1++)
    for(int d2 = 0; d2 < D2; d2++) {
      int ind = d1*D2 +
                d2;
      ptr5D[ind] = ptr4D + ind*D3;
    }
  for(int d1 = 0; d1 < D1; d1++) {
    ptr6D[d1] = ptr5D + d1*D2;
  }
#ifdef _OPENACC
  T * d_ptr = (T*) GetDevicePtr(ptr);
  T ** d_ptr2D = (T**) GetDevicePtr(ptr2D);
  T *** d_ptr3D = (T***) GetDevicePtr(ptr3D);
  T **** d_ptr4D = (T****) GetDevicePtr(ptr4D);
  T ***** d_ptr5D = (T*****) GetDevicePtr(ptr5D);
  T ****** d_ptr6D = (T******) GetDevicePtr(ptr6D);
  T ** tmp2D = (T**) malloc(D1*D2*D3*D4*D5*sizeof(T*));
  T *** tmp3D = (T***) malloc(D1*D2*D3*D4*sizeof(T**));
  T **** tmp4D = (T****) malloc(D1*D2*D3*sizeof(T***));
  T ***** tmp5D = (T*****) malloc(D1*D2*sizeof(T****));
  T ****** tmp6D = (T******) malloc(D1*sizeof(T*****));
  for(int d1 = 0; d1 < D1; d1++)
    for(int d2 = 0; d2 < D2; d2++)
      for(int d3 = 0; d3 < D3; d3++)
        for(int d4 = 0; d4 < D4; d4++)
          for(int d5 = 0; d5 < D5; d5++) {
            int ind = d1*D2*D3*D4*D5 +
                      d2*   D3*D4*D5 +
                      d3*      D4*D5 +
                      d4*         D5 +
                      d5;
            tmp2D[ind] = d_ptr + ind*D6;
          }
  for(int d1 = 0; d1 < D1; d1++)
    for(int d2 = 0; d2 < D2; d2++)
      for(int d3 = 0; d3 < D3; d3++)
        for(int d4 = 0; d4 < D4; d4++) {
          int ind = d1*D2*D3*D4 +
                    d2*   D3*D4 +
                    d3*      D4 +
                    d4;
          tmp3D[ind] = d_ptr2D + ind*D5;
        }
  for(int d1 = 0; d1 < D1; d1++)
    for(int d2 = 0; d2 < D2; d2++)
      for(int d3 = 0; d3 < D3; d3++) {
        int ind = d1*D2*D3 +
                  d2*   D3 +
                  d3;
        tmp4D[ind] = d_ptr3D + ind*D4;
      }
  for(int d1 = 0; d1 < D1; d1++)
    for(int d2 = 0; d2 < D2; d2++) {
      int ind = d1*D2 +
                d2;
      tmp5D[ind] = d_ptr4D + ind*D3;
    }
  for(int d1 = 0; d1 < D1; d1++) {
    tmp6D[d1] = d_ptr5D + d1*D2;
  }
  MemcpyCPUtoGPU(tmp2D, d_ptr2D, D1*D2*D3*D4*D5*sizeof(T*));
  MemcpyCPUtoGPU(tmp3D, d_ptr3D, D1*D2*D3*D4*sizeof(T**));
  MemcpyCPUtoGPU(tmp4D, d_ptr4D, D1*D2*D3*sizeof(T***));
  MemcpyCPUtoGPU(tmp5D, d_ptr5D, D1*D2*sizeof(T****));
  MemcpyCPUtoGPU(tmp6D, d_ptr6D, D1*sizeof(T*****));
  free(tmp2D);
  free(tmp3D);
  free(tmp4D);
  free(tmp5D);
  free(tmp6D);
#endif
  return ptr6D;
} 

template <typename T>
void Free6D(
  T ****** ptr, std::size_t D1, std::size_t D2, std::size_t D3,
  std::size_t D4, std::size_t D5, std::size_t D6
)
{
  Free(ptr[0][0][0][0][0], D1*D2*D3*D4*D5*D6*sizeof(T));
  Free(ptr[0][0][0][0], D1*D2*D3*D4*D5*sizeof(T*));
  Free(ptr[0][0][0], D1*D2*D3*D4*sizeof(T**));
  Free(ptr[0][0], D1*D2*D3*sizeof(T***));
  Free(ptr[0], D1*D2*sizeof(T****));
  Free(ptr, D1*sizeof(T*****));
}

template <typename T>
void UpdateGPU6D(
  T ****** ptr, std::size_t D1, std::size_t D2, std::size_t D3,
  std::size_t D4, std::size_t D5, std::size_t D6
)
{
  UpdateGPU(ptr[0][0][0][0][0], D1*D2*D3*D4*D5*D6*sizeof(T));
}

template <typename T>
void UpdateCPU6D(
  T ****** ptr, std::size_t D1, std::size_t D2, std::size_t D3,
  std::size_t D4, std::size_t D5, std::size_t D6
)
{
  UpdateCPU(ptr[0][0][0][0][0], D1*D2*D3*D4*D5*D6*sizeof(T));
}

//////////////////////////////////////////////////////////////////
// T I M I N G ///////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

void StartTimer();
void StopTimer(std::string);
void PrintTimingInformation();

} // end namespace ACCH

#define PGI_COMPARE(ptr,type,size,name,file,func,line)

#if defined(PGICOMPARE) && defined(__PGI)
#define COMPARE(ptr,type,size,name,file,func,line) \
	ACCH::Compare(ptr,type,size,name,file,func,line);
#define GPU_COMPARE(ptr,type,size,numbytes,name,file,func,line) \
	ACCH::GPUCompare(ptr,type,size,numbytes,name,file,func,line);
#else
#define COMPARE(ptr,type,size,name,file,func,line)
#define GPU_COMPARE(ptr,type,size,numbytes,name,file,func,line)
#endif

#if defined(TIMING)
#define START_TIMER ACCH::StartTimer();
#define STOP_TIMER(func) \
	ACCH::StopTimer(func);
#define PRINT_TIMER ACCH::PrintTimingInformation();
#else
#define START_TIMER
#define STOP_TIMER
#define PRINT_TIMER
#endif

#endif






