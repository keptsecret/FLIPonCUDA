#ifndef FOC_GPU_UTIL_H
#define FOC_GPU_UTIL_H

#include "foc.h"

#include <cuda.h>
#include <cuda_runtime_api.h>

#define CUDA_CHECK(EXPR) \
	{ cudaCheck((EXPR), __FILE__, __LINE__); }

inline void cudaCheck(cudaError_t expr, const char* file, int line) {
	if (expr != cudaSuccess) {
		cudaError_t error = cudaGetLastError();
		fprintf(stderr, "CUDA error: %s %s %d\n", cudaGetErrorString(error), file, line);
	}
}

namespace cuda {

inline size_t getCUDADeviceGlobalMem() {
	cudaDeviceProp prop{};
	CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
	return prop.totalGlobalMem;
}

} // namespace cuda

#endif // FOC_GPU_UTIL_H
