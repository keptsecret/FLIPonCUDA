#include "kernels.h"

namespace cuda {

struct KernelCoefficients {
	float coef1, coef2, coef3;
};

FOC_GPU void load_local_memory(int num_values, float* global_data, float* shared_data) {
	size_t tid = threadIdx.x;
	size_t bid = blockIdx.x;

	int block_size = blockDim.x;
	int num_read = ceil((float)num_values / (float)block_size);

	int max_read_thread_id = floor((float)num_values / (float)num_read) - 1;
	int num_remainder = num_values - (max_read_thread_id + 1) * num_read;

	if (tid <= max_read_thread_id) {
		int thread_offset = num_read * tid;
		int global_offset = bid * num_values + thread_offset;

		for (int i = 0; i < num_read; i++) {
			shared_data[thread_offset + i] = global_data[global_offset + i];
		}

		if (tid == max_read_thread_id) {
			thread_offset = thread_offset + num_read;
			global_offset = global_offset + num_read;
			for (int j = 0; j < num_remainder; j++) {
				shared_data[thread_offset + j] = global_data[global_offset + j];
			}
		}
	}

	__syncthreads();
}

FOC_GPU int3 flat_to_3d_index(int flatidx, int isize, int jsize) {
	int i = flatidx % isize;
	int j = (flatidx / isize) % jsize;
	int k = flatidx / (jsize * isize);
	return make_int3(i, j, k);
}

FOC_GPU KernelCoefficients calculate_kernel_coefficients(float r) {
	KernelCoefficients coefs;
	coefs.coef1 = (4.0f / 9.0f) * (1.0f / (r * r * r * r * r * r));
	coefs.coef2 = (17.0f / 9.0f) * (1.0f / (r * r * r * r));
	coefs.coef3 = (22.0f / 9.0f) * (1.0f / (r * r));

	return coefs;
}

FOC_GPU float evaluate_kernel(float rsq, struct KernelCoefficients* coefs) {
	return 1.0 - (*coefs).coef1 * rsq * rsq * rsq + (*coefs).coef2 * rsq * rsq - (*coefs).coef3 * rsq;
}

__global__ void compute_scalar_field_point_value_kernel(float* point_values,
		float* field_data,
		int* batch_offsets,
		int* num_points,
		float* radius,
		float* dx) {
	extern __shared__ float local_points[];

	unsigned int tid = threadIdx.x;
	unsigned int bid = blockIdx.x;
	unsigned int gid = tid + bid * blockDim.x;

	int block_size = blockDim.x;
	int batch_dim = floor(cbrt((float)block_size));
	int num_cells = batch_dim * batch_dim * batch_dim;

	int shared_data_size = 4 * *num_points;
	load_local_memory(shared_data_size, point_values, local_points);

	if (tid >= num_cells) {
		return;
	}

	int3 cell_index = flat_to_3d_index(tid, batch_dim, batch_dim);
	float3 cell_center = make_float3(((float)cell_index.x + 0.5f) * *dx,
			((float)cell_index.y + 0.5f) * *dx,
			((float)cell_index.z + 0.5f) * *dx);

	int3 batch_index = make_int3(batch_offsets[3 * bid + 0],
			batch_offsets[3 * bid + 1],
			batch_offsets[3 * bid + 2]);

	float3 position_offset = make_float3(batch_index.x * batch_dim * *dx,
			batch_index.y * batch_dim * *dx,
			batch_index.z * batch_dim * *dx);

	KernelCoefficients coefficients = calculate_kernel_coefficients(*radius);

	float maxsq = *radius * *radius;
	float sum = 0.0f;
	for (int i = 0; i < shared_data_size; i += 4) {
		float3 p = make_float3(local_points[i + 0] - position_offset.x,
				local_points[i + 1] - position_offset.y,
				local_points[i + 2] - position_offset.z);
		float value = local_points[i + 3];

		float3 vector = make_float3(p.x - cell_center.x,
				p.y - cell_center.y,
				p.z - cell_center.z);
		float rsq = vector.x * vector.x + vector.y * vector.y + vector.z * vector.z;

		if (rsq < maxsq) {
			sum += value * evaluate_kernel(rsq, &coefficients);
		}
	}

	int fieldidx = bid * num_cells + tid;
	field_data[fieldidx] = sum;
}

__global__ void compute_scalar_weight_field_point_value_kernel(float* point_values,
		float* field_data,
		int* batch_offsets,
		int* num_points,
		int* num_batches,
		float* radius,
		float* dx) {
	extern __shared__ float local_points[];

	unsigned int tid = threadIdx.x;
	unsigned int bid = blockIdx.x;
	unsigned int gid = tid + bid * blockDim.x;

	int block_size = blockDim.x;
	int batch_dim = floor(cbrt((float)block_size));
	int num_cells = batch_dim * batch_dim * batch_dim;

	int shared_data_size = 4 * *num_points;
	load_local_memory(shared_data_size, point_values, local_points);

	if (tid >= num_cells) {
		return;
	}

	int3 cell_index = flat_to_3d_index(tid, batch_dim, batch_dim);
	float3 cell_center = make_float3(((float)cell_index.x + 0.5f) * *dx,
			((float)cell_index.y + 0.5f) * *dx,
			((float)cell_index.z + 0.5f) * *dx);

	int3 batch_index = make_int3(batch_offsets[3 * bid + 0],
			batch_offsets[3 * bid + 1],
			batch_offsets[3 * bid + 2]);

	float3 position_offset = make_float3(batch_index.x * batch_dim * *dx,
			batch_index.y * batch_dim * *dx,
			batch_index.z * batch_dim * *dx);

	KernelCoefficients coefficients = calculate_kernel_coefficients(*radius);

	float maxsq = *radius * *radius;
	float fieldsum = 0.0f;
	float weightsum = 0.0f;
	for (int i = 0; i < shared_data_size; i += 4) {
		float3 p = make_float3(local_points[i + 0] - position_offset.x,
				local_points[i + 1] - position_offset.y,
				local_points[i + 2] - position_offset.z);
		float value = local_points[i + 3];

		float3 vector = make_float3(p.x - cell_center.x,
				p.y - cell_center.y,
				p.z - cell_center.z);
		float rsq = vector.x * vector.x + vector.y * vector.y + vector.z * vector.z;

		if (rsq < maxsq) {
			float weight = evaluate_kernel(rsq, &coefficients);
			fieldsum += value * weight;
			weightsum += weight;
		}
	}

	int fieldidx = bid * num_cells + tid;
	int weightfieldidx = *num_batches * num_cells + fieldidx;
	field_data[fieldidx] = fieldsum;
	field_data[weightfieldidx] = weightsum;
}

void launch_scalar_field_point_value_kernel(dim3 blockSize, dim3 gridSize,
		std::vector<float>& pointvalues, std::vector<float>& vfieldh, std::vector<foc::Point3i>& batch_offsets,
		int num_points, float radius, float cellsize) {

	// copy data to device
	float *pointvalues_d, *vfield_d;
	int* offsets_d;
	int *num_points_d, *num_batches_d;
	float *radius_d, *cellsize_d;

	CUDA_CHECK(cudaMalloc(&pointvalues_d, sizeof(float) * pointvalues.size()));
	CUDA_CHECK(cudaMalloc(&vfield_d, sizeof(float) * vfieldh.size()));
	CUDA_CHECK(cudaMalloc(&offsets_d, 3 * sizeof(int) * batch_offsets.size()));
	CUDA_CHECK(cudaMalloc(&num_points_d, sizeof(int)));
	CUDA_CHECK(cudaMalloc(&radius_d, sizeof(float)));
	CUDA_CHECK(cudaMalloc(&cellsize_d, sizeof(float)));

	CUDA_CHECK(cudaMemcpy(pointvalues_d, pointvalues.data(), sizeof(float) * pointvalues.size(), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(vfield_d, vfieldh.data(), sizeof(float) * vfieldh.size(), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(offsets_d, batch_offsets.data(), 3 * sizeof(int) * batch_offsets.size(), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(num_points_d, &num_points, sizeof(int), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(radius_d, &radius, sizeof(float), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(cellsize_d, &cellsize, sizeof(float), cudaMemcpyHostToDevice));

	// launch kernel
	int shared_mem_size = num_points * 4 * sizeof(float);
	compute_scalar_field_point_value_kernel<<<gridSize, blockSize, shared_mem_size>>>(pointvalues_d, vfield_d, offsets_d, num_points_d, radius_d, cellsize_d);
	CUDA_CHECK(cudaPeekAtLastError());
	CUDA_CHECK(cudaDeviceSynchronize());

	// copy data from device
	CUDA_CHECK(cudaMemcpy(vfieldh.data(), vfield_d, sizeof(float) * vfieldh.size(), cudaMemcpyDeviceToHost));

	// free device memory
	CUDA_CHECK(cudaFree(pointvalues_d));
	CUDA_CHECK(cudaFree(vfield_d));
	CUDA_CHECK(cudaFree(offsets_d));
	CUDA_CHECK(cudaFree(num_points_d));
	CUDA_CHECK(cudaFree(radius_d));
	CUDA_CHECK(cudaFree(cellsize_d));
}

void launch_scalar_weight_field_point_value_kernel(dim3 blockSize, dim3 gridSize,
		std::vector<float>& pointvalues, std::vector<float>& vfieldh, std::vector<foc::Point3i>& batch_offsets,
		int num_points, float radius, float cellsize) {
	int num_batches = batch_offsets.size();

	// copy data to device
	float *pointvalues_d, *vfield_d;
	int* offsets_d;
	int *num_points_d, *num_batches_d;
	float *radius_d, *cellsize_d;

	CUDA_CHECK(cudaMalloc(&pointvalues_d, sizeof(float) * pointvalues.size()));
	CUDA_CHECK(cudaMalloc(&vfield_d, sizeof(float) * vfieldh.size()));
	CUDA_CHECK(cudaMalloc(&offsets_d, 3 * sizeof(int) * batch_offsets.size()));
	CUDA_CHECK(cudaMalloc(&num_points_d, sizeof(int)));
	CUDA_CHECK(cudaMalloc(&num_batches_d, sizeof(int)));
	CUDA_CHECK(cudaMalloc(&radius_d, sizeof(float)));
	CUDA_CHECK(cudaMalloc(&cellsize_d, sizeof(float)));

	CUDA_CHECK(cudaMemcpy(pointvalues_d, pointvalues.data(), sizeof(float) * pointvalues.size(), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(vfield_d, vfieldh.data(), sizeof(float) * vfieldh.size(), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(offsets_d, batch_offsets.data(), 3 * sizeof(int) * batch_offsets.size(), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(num_points_d, &num_points, sizeof(int), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(num_batches_d, &num_batches, sizeof(int), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(radius_d, &radius, sizeof(float), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(cellsize_d, &cellsize, sizeof(float), cudaMemcpyHostToDevice));

	// launch kernel
	int shared_mem_size = num_points * 4 * sizeof(float);
	compute_scalar_weight_field_point_value_kernel<<<gridSize, blockSize, shared_mem_size>>>(pointvalues_d, vfield_d, offsets_d, num_points_d, num_batches_d, radius_d, cellsize_d);
	CUDA_CHECK(cudaPeekAtLastError());
	CUDA_CHECK(cudaDeviceSynchronize());

	// copy data from device
	CUDA_CHECK(cudaMemcpy(vfieldh.data(), vfield_d, sizeof(float) * vfieldh.size(), cudaMemcpyDeviceToHost));

	// free device memory
	CUDA_CHECK(cudaFree(pointvalues_d));
	CUDA_CHECK(cudaFree(vfield_d));
	CUDA_CHECK(cudaFree(offsets_d));
	CUDA_CHECK(cudaFree(num_points_d));
	CUDA_CHECK(cudaFree(num_batches_d));
	CUDA_CHECK(cudaFree(radius_d));
	CUDA_CHECK(cudaFree(cellsize_d));
}

} // namespace cuda