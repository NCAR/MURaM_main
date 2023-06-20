#ifndef __MURAMACC__
#define __MURAMACC__

#define STRINGIFY(a) #a

// NVIDIA NV TOOLS Profiler Push/Pop
#ifdef USE_NVTX
#include "nvToolsExt.h"
const uint32_t ProfilerColors[] = { 0xff00ff00, 0xff0000ff, 0xffffff00, 0xffff00ff, 0xff00ffff, 0xffff0000, 0xffffffff };
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

// PGI acc_compare pgi_compare for C++
#if defined(DEBUG) && defined(PGICOMPARE) && defined(__PGI)
#include <openacc.h>
#define PGI_COMPARE(ptr,type,numele,name,file,func,line) {          \
    type* pgi_compare_ptr = (type*)ptr;                             \
    if( acc_is_present( pgi_compare_ptr , numele*sizeof(type) ) ) { \
    _Pragma(STRINGIFY(acc update self(pgi_compare_ptr[:numele])))   \
    }                                                               \
    pgi_compare(pgi_compare_ptr, #type, numele,                     \
    name, file, func, line);                                        \
    if( acc_is_present( pgi_compare_ptr , numele*sizeof(type) ) ) { \
    _Pragma(STRINGIFY(acc update device(pgi_compare_ptr[:numele])))   \
    }                                                               \
}
#else
#define PGI_COMPARE(ptr,type,size,name,file,func,line)
#endif

#if defined(DEBUG) && defined(ACCCOMPARE) && defined(_OPENACC) && defined(__PGI)
#define ACC_COMPARE(ptr,type,size) {                      \
    type* acc_compare_ptr = (type*) ptr;                  \
    acc_compare(acc_compare_ptr, size);                   \
}
#else
#define ACC_COMPARE(ptr,type,size)
#endif

#ifdef _OPENACC
#include <openacc.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <cstdint>

void * accMalloc(size_t numbytes);
void accFree(void * ptr);
void * gpuMalloc(size_t numbytes);
void gpuFree(void * ptr);

#endif























