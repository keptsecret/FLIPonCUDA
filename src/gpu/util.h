#ifndef FOC_GPU_UTIL_H
#define FOC_GPU_UTIL_H

#include "foc.h"

#include <cuda.h>
#include <cuda_runtime_api.h>

#define CUDA_CHECK(EXPR)                                        \
    if (EXPR != cudaSuccess) {                                  \
        cudaError_t error = cudaGetLastError();                 \
        printf("CUDA error: %s", cudaGetErrorString(error)); 	\
    } else /* eat semicolon */

namespace cuda {

inline size_t getCUDADeviceGlobalMem() {
	cudaDeviceProp prop{};
	CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
	return prop.totalGlobalMem;
}

}

#endif // FOC_GPU_UTIL_H
