#include "kernels.h"

#define BATCH_DIM 5

#define U_OFFSET 0
// V_OFFSET + (CHUNK_WIDTH+4)*(CHUNK_HEIGHT+3)*(CHUNK_DEPTH+4)
#define V_OFFSET 648
// W_OFFSET + (CHUNK_WIDTH+4)*(CHUNK_HEIGHT+4)*(CHUNK_DEPTH+3)
#define W_OFFSET 1296
#define VFIELD_SIZE 1944
// VFIELD_SIZE / 4
#define MAX_VFIELD_LOAD_LOCAL_ID 486

namespace cuda {

FOC_GPU float cubic_interpolate(float p[4], float x) {
	return p[1] + 0.5 * x * (p[2] - p[0] + x * (2.0 * p[0] - 5.0 * p[1] + 4.0 * p[2] - p[3] + x * (3.0 * (p[1] - p[2]) + p[3] - p[0])));
}

FOC_GPU float bicubic_interpolate(float p[4][4], float x, float y) {
	float arr[4];
	arr[0] = cubic_interpolate(p[0], x);
	arr[1] = cubic_interpolate(p[1], x);
	arr[2] = cubic_interpolate(p[2], x);
	arr[3] = cubic_interpolate(p[3], x);
	return cubic_interpolate(arr, y);
}

FOC_GPU float tricubic_interpolate(float p[4][4][4], float x, float y, float z) {
	float arr[4];
	arr[0] = bicubic_interpolate(p[0], x, y);
	arr[1] = bicubic_interpolate(p[1], x, y);
	arr[2] = bicubic_interpolate(p[2], x, y);
	arr[3] = bicubic_interpolate(p[3], x, y);
	return cubic_interpolate(arr, z);
}

FOC_GPU int flatten_index(int i, int j, int k, int isize, int jsize) {
	return i + isize * (j + jsize * k);
}

FOC_GPU void fill_interpolation_data(float* vfield,
		int3 voffset, int vwidth, int vheight,
		float points[4][4][4]) {
	int flatidx;
	flatidx = flatten_index(voffset.x, voffset.y + 0, voffset.z + 0, vwidth, vheight);
	points[0][0][0] = vfield[flatidx + 0];
	points[0][0][1] = vfield[flatidx + 1];
	points[0][0][2] = vfield[flatidx + 2];
	points[0][0][3] = vfield[flatidx + 3];

	flatidx = flatten_index(voffset.x, voffset.y + 1, voffset.z + 0, vwidth, vheight);
	points[0][1][0] = vfield[flatidx + 0];
	points[0][1][1] = vfield[flatidx + 1];
	points[0][1][2] = vfield[flatidx + 2];
	points[0][1][3] = vfield[flatidx + 3];

	flatidx = flatten_index(voffset.x, voffset.y + 2, voffset.z + 0, vwidth, vheight);
	points[0][2][0] = vfield[flatidx + 0];
	points[0][2][1] = vfield[flatidx + 1];
	points[0][2][2] = vfield[flatidx + 2];
	points[0][2][3] = vfield[flatidx + 3];

	flatidx = flatten_index(voffset.x, voffset.y + 3, voffset.z + 0, vwidth, vheight);
	points[0][3][0] = vfield[flatidx + 0];
	points[0][3][1] = vfield[flatidx + 1];
	points[0][3][2] = vfield[flatidx + 2];
	points[0][3][3] = vfield[flatidx + 3];

	flatidx = flatten_index(voffset.x, voffset.y + 0, voffset.z + 1, vwidth, vheight);
	points[1][0][0] = vfield[flatidx + 0];
	points[1][0][1] = vfield[flatidx + 1];
	points[1][0][2] = vfield[flatidx + 2];
	points[1][0][3] = vfield[flatidx + 3];

	flatidx = flatten_index(voffset.x, voffset.y + 1, voffset.z + 1, vwidth, vheight);
	points[1][1][0] = vfield[flatidx + 0];
	points[1][1][1] = vfield[flatidx + 1];
	points[1][1][2] = vfield[flatidx + 2];
	points[1][1][3] = vfield[flatidx + 3];

	flatidx = flatten_index(voffset.x, voffset.y + 2, voffset.z + 1, vwidth, vheight);
	points[1][2][0] = vfield[flatidx + 0];
	points[1][2][1] = vfield[flatidx + 1];
	points[1][2][2] = vfield[flatidx + 2];
	points[1][2][3] = vfield[flatidx + 3];

	flatidx = flatten_index(voffset.x, voffset.y + 3, voffset.z + 1, vwidth, vheight);
	points[1][3][0] = vfield[flatidx + 0];
	points[1][3][1] = vfield[flatidx + 1];
	points[1][3][2] = vfield[flatidx + 2];
	points[1][3][3] = vfield[flatidx + 3];

	flatidx = flatten_index(voffset.x, voffset.y + 0, voffset.z + 2, vwidth, vheight);
	points[2][0][0] = vfield[flatidx + 0];
	points[2][0][1] = vfield[flatidx + 1];
	points[2][0][2] = vfield[flatidx + 2];
	points[2][0][3] = vfield[flatidx + 3];

	flatidx = flatten_index(voffset.x, voffset.y + 1, voffset.z + 2, vwidth, vheight);
	points[2][1][0] = vfield[flatidx + 0];
	points[2][1][1] = vfield[flatidx + 1];
	points[2][1][2] = vfield[flatidx + 2];
	points[2][1][3] = vfield[flatidx + 3];

	flatidx = flatten_index(voffset.x, voffset.y + 2, voffset.z + 2, vwidth, vheight);
	points[2][2][0] = vfield[flatidx + 0];
	points[2][2][1] = vfield[flatidx + 1];
	points[2][2][2] = vfield[flatidx + 2];
	points[2][2][3] = vfield[flatidx + 3];

	flatidx = flatten_index(voffset.x, voffset.y + 3, voffset.z + 2, vwidth, vheight);
	points[2][3][0] = vfield[flatidx + 0];
	points[2][3][1] = vfield[flatidx + 1];
	points[2][3][2] = vfield[flatidx + 2];
	points[2][3][3] = vfield[flatidx + 3];

	flatidx = flatten_index(voffset.x, voffset.y + 0, voffset.z + 3, vwidth, vheight);
	points[3][0][0] = vfield[flatidx + 0];
	points[3][0][1] = vfield[flatidx + 1];
	points[3][0][2] = vfield[flatidx + 2];
	points[3][0][3] = vfield[flatidx + 3];

	flatidx = flatten_index(voffset.x, voffset.y + 1, voffset.z + 3, vwidth, vheight);
	points[3][1][0] = vfield[flatidx + 0];
	points[3][1][1] = vfield[flatidx + 1];
	points[3][1][2] = vfield[flatidx + 2];
	points[3][1][3] = vfield[flatidx + 3];

	flatidx = flatten_index(voffset.x, voffset.y + 2, voffset.z + 3, vwidth, vheight);
	points[3][2][0] = vfield[flatidx + 0];
	points[3][2][1] = vfield[flatidx + 1];
	points[3][2][2] = vfield[flatidx + 2];
	points[3][2][3] = vfield[flatidx + 3];

	flatidx = flatten_index(voffset.x, voffset.y + 3, voffset.z + 3, vwidth, vheight);
	points[3][3][0] = vfield[flatidx + 0];
	points[3][3][1] = vfield[flatidx + 1];
	points[3][3][2] = vfield[flatidx + 2];
	points[3][3][3] = vfield[flatidx + 3];
}

FOC_GPU float interpolate_U(float3 pos, float dx, float invdx, float* ufield) {
	pos.y -= 0.5 * dx;
	pos.z -= 0.5 * dx;

	int3 index = make_int3(floor(pos.x * invdx),
			floor(pos.y * invdx),
			floor(pos.z * invdx));

	float3 index_offset = make_float3(index.x * dx,
			index.y * dx,
			index.z * dx);

	float3 interp_pos = make_float3(invdx * (pos.x - index_offset.x),
			invdx * (pos.y - index_offset.y),
			invdx * (pos.z - index_offset.z));

	int3 vfield_index_offset = make_int3(index.x - 1 + 1,
			index.y - 1 + 2,
			index.z - 1 + 2);

	float points[4][4][4];
	int vwidth = BATCH_DIM + 3;
	int vheight = BATCH_DIM + 4;

	fill_interpolation_data(ufield, vfield_index_offset, vwidth, vheight, points);

	return tricubic_interpolate(points, interp_pos.x,
			interp_pos.y,
			interp_pos.z);
}

FOC_GPU float interpolate_V(float3 pos, float dx, float invdx, float* vfield) {
	pos.x -= 0.5 * dx;
	pos.z -= 0.5 * dx;

	int3 index = make_int3(floor(pos.x * invdx),
			floor(pos.y * invdx),
			floor(pos.z * invdx));

	float3 index_offset = make_float3(index.x * dx,
			index.y * dx,
			index.z * dx);

	float3 interp_pos = make_float3(invdx * (pos.x - index_offset.x),
			invdx * (pos.y - index_offset.y),
			invdx * (pos.z - index_offset.z));

	int3 vfield_index_offset = make_int3(index.x - 1 + 2,
			index.y - 1 + 1,
			index.z - 1 + 2);

	float points[4][4][4];
	int vwidth = BATCH_DIM + 4;
	int vheight = BATCH_DIM + 3;

	fill_interpolation_data(vfield, vfield_index_offset, vwidth, vheight, points);

	return tricubic_interpolate(points, interp_pos.x,
			interp_pos.y,
			interp_pos.z);
}

FOC_GPU float interpolate_W(float3 pos, float dx, float invdx, float* wfield) {
	pos.x -= 0.5 * dx;
	pos.y -= 0.5 * dx;

	int3 index = make_int3(floor(pos.x * invdx),
			floor(pos.y * invdx),
			floor(pos.z * invdx));

	float3 index_offset = make_float3(index.x * dx,
			index.y * dx,
			index.z * dx);

	float3 interp_pos = make_float3(invdx * (pos.x - index_offset.x),
			invdx * (pos.y - index_offset.y),
			invdx * (pos.z - index_offset.z));

	int3 vfield_index_offset = make_int3(index.x - 1 + 2,
			index.y - 1 + 2,
			index.z - 1 + 1);

	float points[4][4][4];
	int vwidth = BATCH_DIM + 4;
	int vheight = BATCH_DIM + 4;

	fill_interpolation_data(wfield, vfield_index_offset, vwidth, vheight, points);

	return tricubic_interpolate(points, interp_pos.x,
			interp_pos.y,
			interp_pos.z);
}

__global__ void tricubic_interpolate_kernel(float* particles, float* vfielddata, int* batch_offsets, float* cellsize) {
	// This kernel essentially implements MACGrid::getVelocityAt() on GPU

	__shared__ float vfield[VFIELD_SIZE];

	unsigned int tid = threadIdx.x;
	unsigned int bid = blockIdx.x;
	unsigned int gid = tid + bid * blockDim.x;

	if (tid < MAX_VFIELD_LOAD_LOCAL_ID) {
		int local_offset = 4 * tid;
		int vfield_data_offset = bid * VFIELD_SIZE + local_offset;

		vfield[local_offset + 0] = vfielddata[vfield_data_offset + 0];
		vfield[local_offset + 1] = vfielddata[vfield_data_offset + 1];
		vfield[local_offset + 2] = vfielddata[vfield_data_offset + 2];
		vfield[local_offset + 3] = vfielddata[vfield_data_offset + 3];
	}

	__syncthreads();

	float3 pos = make_float3(particles[3 * gid + 0],
			particles[3 * gid + 1],
			particles[3 * gid + 2]);

	int3 batch_offset = make_int3(batch_offsets[3 * bid + 0],
			batch_offsets[3 * bid + 1],
			batch_offsets[3 * bid + 2]);

	int3 index_offset = make_int3(batch_offset.x * BATCH_DIM,
			batch_offset.y * BATCH_DIM,
			batch_offset.z * BATCH_DIM);

	float3 pos_offset = make_float3((float)index_offset.x * *cellsize,
			(float)index_offset.y * *cellsize,
			(float)index_offset.z * *cellsize);

	float3 local_pos = make_float3(pos.x - pos_offset.x,
			pos.y - pos_offset.y,
			pos.z - pos_offset.z);

	float invcs = 1.0 / *cellsize;
	float x = interpolate_U(local_pos, *cellsize, invcs, &(vfield[U_OFFSET]));
	float y = interpolate_V(local_pos, *cellsize, invcs, &(vfield[V_OFFSET]));
	float z = interpolate_W(local_pos, *cellsize, invcs, &(vfield[W_OFFSET]));

	particles[3 * gid + 0] = x;
	particles[3 * gid + 1] = y;
	particles[3 * gid + 2] = z;
}

void launch_tricubicinterpolate_kernel(dim3 blockSize, dim3 gridSize,
		std::vector<foc::Vector3f>& particles, std::vector<float>& vfield, std::vector<foc::Point3i>& batch_offsets, float cellsize,
		std::vector<foc::Vector3f>& outputs) {
	// copy data to device
	float *particles_d, *vfield_d;
	int* offsets_d;
	float* cellsize_d;

	CUDA_CHECK(cudaMalloc(&particles_d, 3 * sizeof(float) * particles.size()));
	CUDA_CHECK(cudaMalloc(&vfield_d, sizeof(float) * vfield.size()));
	CUDA_CHECK(cudaMalloc(&offsets_d, 3 * sizeof(int) * batch_offsets.size()));
	CUDA_CHECK(cudaMalloc(&cellsize_d, sizeof(float)));

	CUDA_CHECK(cudaMemcpy(particles_d, particles.data(), 3 * sizeof(float) * particles.size(), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(vfield_d, vfield.data(), sizeof(float) * vfield.size(), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(offsets_d, batch_offsets.data(), 3 * sizeof(int) * batch_offsets.size(), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(cellsize_d, &cellsize, sizeof(float), cudaMemcpyHostToDevice));

	// launch kernel
	tricubic_interpolate_kernel<<<gridSize, blockSize>>>(particles_d, vfield_d, offsets_d, cellsize_d);
	CUDA_CHECK(cudaPeekAtLastError());
	CUDA_CHECK(cudaDeviceSynchronize());

	// copy data from device
	outputs.resize(particles.size());
	CUDA_CHECK(cudaMemcpy(outputs.data(), particles_d, 3 * sizeof(float) * particles.size(), cudaMemcpyDeviceToHost));

	// free device memory
	CUDA_CHECK(cudaFree(particles_d));
	CUDA_CHECK(cudaFree(vfield_d));
	CUDA_CHECK(cudaFree(offsets_d));
	CUDA_CHECK(cudaFree(cellsize_d));
}

} // namespace cuda
