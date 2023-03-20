#ifndef FOC_H
#define FOC_H

#if defined(__CUDA_ARCH__)
#define IS_GPU_CODE
#endif

#include <stdint.h>

#if defined(CUDA_ENABLED) && defined(__CUDACC__)
#include <cuda.h>
#include <cuda_runtime_api.h>
#define FOC_CPU_GPU __host__ __device__
#define FOC_GPU __device__
#if defined(IS_GPU_CODE)
#define FOC_CONST __device__ const
#else
#define FOC_CONST const
#endif
#else
#define FOC_CONST const
#define FOC_CPU_GPU
#define FOC_GPU
#endif

#include <float.h>
#include <limits.h>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <new>
#include <string>
#include <thread>
#include <vector>

namespace foc {

// Global forward declarations
template <typename T>
class Vector2;

template <typename T>
class Vector3;

template <typename T>
class Point3;

template <typename T>
class Point2;

using Point2f = Point2<float>;
using Point2i = Point2<int>;
using Point3f = Point3<float>;
using Vector2f = Vector2<float>;
using Vector2i = Vector2<int>;
using Vector3f = Vector3<float>;

template <typename T>
class Array3D;

} // namespace foc

#endif // FOC_H
