#ifndef FOC_GPU_TRICUBICINTERPOLATE_H
#define FOC_GPU_TRICUBICINTERPOLATE_H

#include "foc.h"
#include "util.h"
#include "vecmath.h"

namespace cuda {

extern void launch_tricubicinterpolate_kernel(dim3 blockSize, dim3 gridSize,
		std::vector<foc::Vector3f>& particles, std::vector<float>& vfield, std::vector<foc::Point3i>& batch_offsets, float cellsize,
		std::vector<foc::Vector3f>& output);

extern void launch_scalar_field_point_value_kernel(dim3 blockSize, dim3 gridSize,
		std::vector<float>& pointvalues, std::vector<float>& vfield, std::vector<foc::Point3i>& batch_offsets,
		int num_points, float radius, float cellsize);

extern void launch_scalar_weight_field_point_value_kernel(dim3 blockSize, dim3 gridSize,
		std::vector<float>& pointvalues, std::vector<float>& vfield, std::vector<foc::Point3i>& batch_offsets,
		int num_points, float radius, float cellsize);

} // namespace cuda

#endif // FOC_GPU_TRICUBICINTERPOLATE_H
