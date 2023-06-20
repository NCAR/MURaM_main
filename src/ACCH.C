#include <stdio.h>
#include <cstring>

#include "ACCH.h"

namespace ACCH
{
#ifndef CACHE_SIZE
  int shared_size;
#endif
  int gpu_num;
}

#ifdef ACCDEBUG
double CPUMB = 0;
double GPUMB = 0;
#endif

#if defined(_OPENACC) || defined(__PGI)
#include <openacc.h>
#endif

void * ACCH::Malloc(
  std::size_t numbytes
)
{
  void * ptr = malloc(numbytes);
#ifdef ACCDEBUG
  double MB = (double) numbytes/(1024*1024);
  CPUMB += MB;
  printf("%.1f MB allocated on CPU (total %.1f): %p\n", MB, CPUMB, ptr);
#endif
#ifdef _OPENACC
  void * d_ptr = acc_create(ptr, numbytes);
#ifdef ACCDEBUG
  GPUMB += MB;
  printf("%.1f MB allocated on GPU (total %.1f): %p\n", MB, GPUMB, d_ptr);
#endif
#endif
  return ptr;
}

void ACCH::Free(
  void * ptr, std::size_t numbytes
)
{
#ifdef ACCDEBUG
  double MB = (double) numbytes/(1024*1024);
#endif
#ifdef _OPENACC
#ifdef ACCDEBUG
  GPUMB -= MB;
  void * d_ptr = ACCH::GetDevicePtr(ptr);
  printf("%.1f MB deallocated from GPU (total %.1f): %p\n", MB, GPUMB, d_ptr);
#endif
  acc_delete(ptr, numbytes);
#endif
#ifdef ACCDEBUG
  CPUMB -= MB;
  printf("%.1f MB deallocated from CPU (total %.1f): %p\n", MB, CPUMB, ptr);
#endif
  free(ptr);
}

void * ACCH::GPUMalloc(
  std::size_t numbytes
)
{
#ifdef _OPENACC
  void * ptr = acc_malloc(numbytes);
#ifdef ACCDEBUG
  double MB = (double) numbytes/(1024*1024);
  GPUMB += MB;
  printf("%.1f MB allocated on GPU (total %.1f): %p\n", MB, GPUMB, ptr);
#endif
  return ptr;
#else
void * ptr =  malloc(numbytes);
#ifdef ACCDEBUG
  printf("%zu MB allocated on CPU: %p\n", (double) numbytes/(1024*1024), ptr);
#endif
  return ptr;
#endif
}

void ACCH::GPUFree(
  void * ptr
)
{
#ifdef _OPENACC
  acc_free(ptr);
#else
  free(ptr);
#endif
}

void ACCH::GPUFree(
  void * ptr, std::size_t numbytes
)
{
#ifdef _OPENACC
#ifdef ACCDEBUG
  GPUMB -= (double) numbytes/(1024*1024);
#endif
  acc_free(ptr);
#else
#ifdef ACCDEBUG
  CPUMB -= (double) numbytes/(1024*1024);
#endif
  free(ptr);
#endif
}

void ACCH::UpdateCPU(
  void * ptr, std::size_t numbytes
)
{
#ifdef _OPENACC
  acc_update_self(ptr, numbytes);
#endif
}

void ACCH::UpdateGPU(
  void * ptr, std::size_t numbytes
)
{
#ifdef _OPENACC
  acc_update_device(ptr, numbytes);
#endif
}

void ACCH::UpdateCPU(
  void * ptr, std::size_t start, std::size_t numbytes
)
{
#ifdef _OPENACC
  void * d_ptr = ACCH::GetDevicePtr(ptr);
  acc_memcpy_from_device(
    static_cast<char*>(ptr)+start,
    static_cast<char*>(d_ptr)+start,
    numbytes
  );
#endif
}

void ACCH::UpdateGPU(
  void * ptr, std::size_t start, std::size_t numbytes
)
{
#ifdef _OPENACC
  void * d_ptr = ACCH::GetDevicePtr(ptr);
  acc_memcpy_to_device(
    static_cast<char*>(d_ptr)+start,
    static_cast<char*>(ptr)+start,
    numbytes
  );
#endif
}

void ACCH::MemcpyCPUtoGPU(
  void * h_ptr, void * d_ptr, std::size_t numbytes
)
{
#ifdef _OPENACC
  acc_memcpy_to_device(d_ptr, h_ptr, numbytes);
#endif
}

void ACCH::MemcpyGPUtoCPU(
  void * h_ptr, void * d_ptr, std::size_t numbytes
)
{
#ifdef _OPENACC
  acc_memcpy_from_device(h_ptr, d_ptr, numbytes);
#endif
}

void ACCH::MemcpyCPU(
  void * src, void * dest, std::size_t bytes
)
{
  std::memcpy(dest, src, bytes);
}

void ACCH::MemcpyGPU(
  void * src, void * dest, std::size_t bytes
)
{
#ifdef _OPENACC
  acc_memcpy_device(dest, src, bytes);
#else
  std::memcpy(dest, src, bytes);
#endif
}

void ACCH::Create(
  void * ptr, std::size_t bytes
)
{
#ifdef _OPENACC
  acc_create(ptr, bytes);
#endif
}

void ACCH::Copyin(
  void * ptr, std::size_t bytes
)
{
#ifdef _OPENACC
  acc_copyin(ptr, bytes);
#endif
}

void ACCH::Delete(
  void * ptr, std::size_t bytes
)
{
#ifdef _OPENACC
  acc_delete(ptr, bytes);
#endif
}

void ACCH::Copyout(
  void * ptr, std::size_t bytes
)
{
#ifdef _OPENACC
  acc_copyout(ptr, bytes);
#endif
}

bool ACCH::Present(
  void * ptr, std::size_t bytes
)
{
#ifdef _OPENACC
  return acc_is_present(ptr, bytes);
#else
  return true;
#endif
}

void ACCH::Compare(
  void * data, const char * datatype, std::size_t numelements,
  const char * dataname, const char * filename, 
  const char * funcname, unsigned int linenum
)
{
#if defined(__PGI) && defined(PGICOMPARE)
  pgi_compare(
    data, datatype, numelements,
    dataname, filename, 
    funcname, linenum
  );
#endif
}

void ACCH::GPUCompare(
  void * data, const char * datatype, std::size_t numelements, std::size_t numbytes,
  const char * dataname, const char * filename, 
  const char * funcname, unsigned int linenum
)
{
#if defined(__PGI) && defined(PGICOMPARE)
  UpdateCPU(data, numbytes);
  Compare(
    data, datatype, numelements,
    dataname, filename, funcname, linenum
  );
  UpdateGPU(data, numbytes);
#endif
}

int ACCH::GetNumGPU()
{
#ifdef _OPENACC
  return acc_get_num_devices(acc_device_default);
#else
  return 1;
#endif
}

void ACCH::SetGPU(
  int num
)
{
  ACCH::gpu_num = num;
#ifdef _OPENACC
  acc_set_device_num(num, acc_device_default);
#endif
}

int ACCH::GetGPU()
{
  return ACCH::gpu_num;
}

void * ACCH::GetDevicePtr(
  void * h_ptr
)
{
#ifdef _OPENACC
  void * d_ptr = acc_deviceptr(h_ptr);
  if(!d_ptr) {
    printf("Warning GetDevicePtr: data not found on device.\n");
    return h_ptr;
  } else {
    return d_ptr;
  }
#else
  return h_ptr;
#endif
}

///////////////////////////////////////////////////////////////////////////
// T I M I N G ////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#include <unordered_map>
#include <vector>
#include <stack>
#include <mpi.h>

namespace
{
  std::stack<double> startTimers;
  std::unordered_map<std::string, std::vector<double> > times;
};

void ACCH::StartTimer()
{
  double st = MPI_Wtime();
  startTimers.push(st);
}

void ACCH::StopTimer(std::string kernel)
{
  double et = MPI_Wtime();
  double st = startTimers.top();
  startTimers.pop();
  if(times.find(kernel) == times.end()) {
    times.insert(std::make_pair(kernel, std::vector<double>()));
  }
  times.at(kernel).push_back(et-st);
}

void ACCH::PrintTimingInformation()
{
  std::cout << "Collected timing information ***\n";
  for(auto it = times.begin(); it != times.end(); ++it) {
    std::cout << it->first << ": ";
    std::vector<double> runtimes = it->second;
    int count = runtimes.size();
    double time_taken = 0.0;
    for(int i = 0; i < count; i++) {
      time_taken += runtimes[i];
    }
    std::cout << count << " calls " << time_taken/count << " avg " << time_taken << " total\n";
  }
}
































