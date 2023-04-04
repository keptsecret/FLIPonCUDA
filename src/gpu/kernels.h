#ifndef FOC_GPU_TRICUBICINTERPOLATE_H
#define FOC_GPU_TRICUBICINTERPOLATE_H

#include "foc.h"
#include "util.h"
#include "vecmath.h"

namespace cuda {

extern void launch_tricubicinterpolate_kernel(dim3 blockSize, dim3 gridSize,
		std::vector<foc::Vector3f>& particles, std::vector<float>& vfield, std::vector<foc::Point3i>& batch_offsets, float cellsize,
		std::vector<foc::Vector3f>& output);

} // namespace cuda

#endif // FOC_GPU_TRICUBICINTERPOLATE_H
