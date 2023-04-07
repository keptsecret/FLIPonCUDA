#ifndef FOC_H
#define FOC_H

#if defined(__CUDA_ARCH__)
#define IS_GPU_CODE
#endif

#include <stdint.h>
#include <cstddef>

#if defined(FOC_BUILD_GPU) && defined(__CUDACC__)
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
using Point3i = Point3<int>;
using Vector2f = Vector2<float>;
using Vector2i = Vector2<int>;
using Vector3f = Vector3<float>;

template <typename T>
class Bounds3;
using Bounds3f = Bounds3<float>;
using Bounds3i = Bounds3<int>;

template <typename T>
class Array3D;

// Global constants
#ifdef _MSC_VER
#define MaxFloat std::numeric_limits<float>::max()
#define Infinity std::numeric_limits<float>::infinity()
#else
static FOC_CONST float MaxFloat = std::numeric_limits<float>::max();
static FOC_CONST float Infinity = std::numeric_limits<float>::infinity();
#endif
#ifdef _MSC_VER
#define MachineEpsilon (std::numeric_limits<float>::epsilon() * 0.5)
#else
static constexpr float MachineEpsilon =
		std::numeric_limits<float>::epsilon() * 0.5;
#endif
static FOC_CONST double ShadowEpsilon = 1e-8f;
static FOC_CONST float Pi = 3.14159265358979323846;
static FOC_CONST float InvPi = 0.31830988618379067154;
static FOC_CONST float Inv2Pi = 0.15915494309189533577;
static FOC_CONST float Inv4Pi = 0.07957747154594766788;
static FOC_CONST float PiOver2 = 1.57079632679489661923;
static FOC_CONST float PiOver4 = 0.78539816339744830961;
static FOC_CONST float Sqrt2 = 1.41421356237309504880;

inline double randomDouble(double min_val, double max_val) {
	return min_val + (double)rand() / ((double)RAND_MAX / (max_val - min_val));
}

FOC_CPU_GPU inline float gamma(int n) {
	return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
}

inline double cubicInterpolate(double p[4], double x) {
	return p[1] + 0.5 * x * (p[2] - p[0] + x * (2.0 * p[0] - 5.0 * p[1] + 4.0 * p[2] - p[3] + x * (3.0 * (p[1] - p[2]) + p[3] - p[0])));
}

inline double bicubicInterpolate(double p[4][4], double x, double y) {
	double arr[4];
	arr[0] = cubicInterpolate(p[0], x);
	arr[1] = cubicInterpolate(p[1], x);
	arr[2] = cubicInterpolate(p[2], x);
	arr[3] = cubicInterpolate(p[3], x);
	return cubicInterpolate(arr, y);
}

inline double tricubicInterpolate(double p[4][4][4], double x, double y, double z) {
	double arr[4];
	arr[0] = bicubicInterpolate(p[0], x, y);
	arr[1] = bicubicInterpolate(p[1], x, y);
	arr[2] = bicubicInterpolate(p[2], x, y);
	arr[3] = bicubicInterpolate(p[3], x, y);
	return cubicInterpolate(arr, z);
}

} // namespace foc

#endif // FOC_H
