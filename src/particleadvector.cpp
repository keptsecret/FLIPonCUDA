#include "particleadvector.h"

#include "gpu/kernels.h"

namespace foc {

ParticleAdvector::ParticleAdvector() {
}

void ParticleAdvector::tricubicInterpolate(std::vector<Vector3f>& particles, MACGrid* vfield, std::vector<Vector3f>& output) {
#ifndef FOC_BUILD_GPU
	output.resize(particles.size());

	for (int i = 0; i < particles.size(); i++) {
		output[i] = vfield->getVelocityAt(particles[i].x, particles[i].y, particles[i].z);
	}

	validateOutput(output);
#else
	vfield->getGridDimensions(&isize, &jsize, &ksize);
	cellsize = vfield->getGridCellSize();

	int batchi = batchDim;
	int batchj = batchDim;
	int batchk = batchDim;
	int batchgridi = std::ceil((double)isize / (double)batchi);
	int batchgridj = std::ceil((double)jsize / (double)batchj);
	int batchgridk = std::ceil((double)ksize / (double)batchk);

	Array3D<ParticleBatch> particleGrid(batchgridi, batchgridj, batchgridk);
	getParticleBatchGrid(batchi * cellsize, batchj * cellsize, batchk * cellsize, particles, particleGrid);

	std::vector<BatchData> batchData;
	getBatchData(vfield, particleGrid, batchData);

	int batchSize = getMaxBatchesPerComputation();
	int numComputations = std::ceil((double)batchData.size() / (double)batchSize);

	output.resize(particles.size());

	for (int i = 0; i < numComputations; i++) {
		int beginIdx = i * batchSize;
		int endIdx = beginIdx + batchSize;
		if (endIdx > batchData.size()) {
			endIdx = batchData.size();
		}

		auto begin = batchData.begin() + beginIdx;
		auto end = batchData.begin() + endIdx;

		std::vector<BatchData> batches;
		batches.insert(batches.begin(), begin, end);

		tricubicInterpolateBatch(batches, output);
	}

	validateOutput(output);
#endif
}

void ParticleAdvector::tricubicInterpolate(std::vector<Vector3f>& particles, MACGrid* vfield) {
	tricubicInterpolate(particles, vfield, particles);
}

void ParticleAdvector::advectParticles(std::vector<Point3f>& particles, MACGrid* vfield, double dt, std::vector<Point3f>& output) {
#ifndef FOC_BUILD_GPU
	output.clear();
	output.reserve(particles.size());
	for (int i = 0; i < particles.size(); i++) {
		output.push_back(RK4(particles[i], dt, vfield));
	}
#else
	output.clear();
	output.reserve(particles.size());

	std::vector<Vector3f> temppos;
	temppos.reserve(particles.size());
	for (int i = 0; i < particles.size(); i++) {
		output.push_back(particles[i]);
		temppos.push_back(Vector3f(particles[i]));
	}

	std::vector<Vector3f> tempdata;
	tempdata.reserve(particles.size());

	tricubicInterpolate(temppos, vfield, tempdata);

	float factor = dt / 6.0;
	for (int i = 0; i < tempdata.size(); i++) {
		output[i] += tempdata[i] * factor;
		tempdata[i] = temppos[i] + tempdata[i] * (0.5 * dt);
	}

	tricubicInterpolate(tempdata, vfield);

	for (int i = 0; i < tempdata.size(); i++) {
		output[i] += tempdata[i] * factor * 2.0;
		tempdata[i] = temppos[i] + tempdata[i] * (0.5 * dt);
	}

	tricubicInterpolate(tempdata, vfield);

	for (int i = 0; i < tempdata.size(); i++) {
		output[i] += tempdata[i] * factor * 2.0;
		tempdata[i] = temppos[i] + tempdata[i] * dt;
	}

	tricubicInterpolate(tempdata, vfield);

	for (int i = 0; i < tempdata.size(); i++) {
		output[i] += tempdata[i] * factor;
	}
#endif
}

void ParticleAdvector::validateOutput(std::vector<Vector3f>& output) {
	for (auto& v : output) {
		if (std::isinf(v.x) || std::isnan(v.x) ||
				std::isinf(v.y) || std::isnan(v.y) ||
				std::isinf(v.z) || std::isnan(v.z)) {
			v = Vector3f();
		}
	}
}

Point3f ParticleAdvector::RK4(Point3f p, double dt, MACGrid* vfield) {
	Vector3f k1 = vfield->getVelocityAt(p);
	Vector3f k2 = vfield->getVelocityAt(p + k1 * 0.5 * dt);
	Vector3f k3 = vfield->getVelocityAt(p + k2 * 0.5 * dt);
	Vector3f k4 = vfield->getVelocityAt(p + k3 * dt);

	return p + (k1 + 2.0f * k2 + 2.0f * k3 + k4) * (dt / 6.0f);
}

int ParticleAdvector::getMaxBatchesPerComputation() {
	int batchPositionDataSize = 3 * sizeof(float) * maxThreadsPerBlock;
	int batchVelocityDataSize = 3 * sizeof(float) * (batchDim + 3) * (batchDim + 4) * (batchDim + 4);
	int batchOffsetDataSize = 3 * sizeof(int);
	int totalSize = batchPositionDataSize + batchVelocityDataSize + batchOffsetDataSize;

	size_t maxGlobalMem = cuda::getCUDADeviceGlobalMem();
	return fmin(std::floor((double)maxGlobalMem / (double)totalSize), maxBatchesPerComputation);
}

void ParticleAdvector::getParticleBatchGrid(double bw, double bh, double bd, std::vector<Vector3f>& particles, Array3D<ParticleBatch>& grid) {
	double bboxeps = 0.01 * cellsize;
	Bounds3f bbox(Point3f(bboxeps, bboxeps, bboxeps), Point3f(grid.width * bw - bboxeps, grid.height * bh - bboxeps, grid.depth * bd - bboxeps));

	Array3D<int> count(grid.width, grid.height, grid.depth, 0);
	for (int i = 0; i < particles.size(); i++) {
		Point3f p(particles[i].x, particles[i].y, particles[i].z);

		if (!bInside(p, bbox)) {
			p = bGetNearestPointInsideBounds(p, ShadowEpsilon, bbox);
		}

		int pi = (int)(p.x / bw);
		int pj = (int)(p.y / bh);
		int pk = (int)(p.z / bd);

		count.set(pi, pj, pk, count(pi, pj, pk) + 1);
	}

	for (int k = 0; k < grid.depth; k++) {
		for (int j = 0; j < grid.height; j++) {
			for (int i = 0; i < grid.width; i++) {
				ParticleBatch* pb = grid.getPointer(i, j, k);
				pb->particles.reserve(count(i, j, k));
				pb->refs.reserve(count(i, j, k));
			}
		}
	}

	for (int i = 0; i < particles.size(); i++) {
		Point3f p(particles[i].x, particles[i].y, particles[i].z);

		if (!bInside(p, bbox)) {
			p = bGetNearestPointInsideBounds(p, ShadowEpsilon, bbox);
		}

		int pi = (int)(p.x / bw);
		int pj = (int)(p.y / bh);
		int pk = (int)(p.z / bd);

		ParticleBatch* pb = grid.getPointer(pi, pj, pk);
		pb->particles.push_back(Vector3f(p));
		pb->refs.push_back(i);
	}

	for (int k = 0; k < grid.depth; k++) {
		for (int j = 0; j < grid.height; j++) {
			for (int i = 0; i < grid.width; i++) {
				Point3f p_min(i * bw, j * bh, k * bd);
				bbox = Bounds3f(p_min + Vector3f(bboxeps, bboxeps, bboxeps),
						p_min + Vector3f(bw - bboxeps, bh - bboxeps, bd - bboxeps));

				ParticleBatch* pb = grid.getPointer(i, j, k);
				for (int idx = 0; idx < pb->particles.size(); idx++) {
					Point3f p = Point3f(pb->particles[idx].x, pb->particles[idx].y, pb->particles[idx].z);
					if (!bInside(p, bbox)) {
						p = bGetNearestPointInsideBounds(p, ShadowEpsilon, bbox);
						pb->particles[idx] = Vector3f(p);
					}
				}
			}
		}
	}
}

void ParticleAdvector::getBatchData(MACGrid* vfield, Array3D<ParticleBatch>& particleGrid, std::vector<BatchData>& batchData) {
	for (int k = 0; k < particleGrid.depth; k++) {
		for (int j = 0; j < particleGrid.height; j++) {
			for (int i = 0; i < particleGrid.width; i++) {
				Point3i bindex(i, j, k);
				ParticleBatch* pb = particleGrid.getPointer(i, j, k);

				if (pb->particles.size() == 0) {
					continue;
				}

				Point3i indexOffset(bindex.x * batchDim, bindex.y * batchDim, bindex.z * batchDim);
				double dx = vfield->getGridCellSize();
				Point3f positionOffset = gridIndexToPosition(indexOffset.x, indexOffset.y, indexOffset.z, dx);

				Array3D<float>* ugrid = vfield->getArrayU();
				Array3D<float>* vgrid = vfield->getArrayV();
				Array3D<float>* wgrid = vfield->getArrayW();

				Point3i ugridOffset(indexOffset.x - 1, indexOffset.y - 2, indexOffset.z - 2);
				Point3i vgridOffset(indexOffset.x - 2, indexOffset.y - 1, indexOffset.z - 2);
				Point3i wgridOffset(indexOffset.x - 2, indexOffset.y - 2, indexOffset.z - 1);

				Array3DView<float> ugridView(batchDim + 3, batchDim + 4, batchDim + 4, ugridOffset.x, ugridOffset.y, ugridOffset.z, ugrid);
				Array3DView<float> vgridView(batchDim + 4, batchDim + 3, batchDim + 4, vgridOffset.x, vgridOffset.y, vgridOffset.z, vgrid);
				Array3DView<float> wgridView(batchDim + 4, batchDim + 4, batchDim + 3, wgridOffset.x, wgridOffset.y, wgridOffset.z, wgrid);

				int blockSize = maxThreadsPerBlock;
				int numBatches = std::ceil((double)pb->particles.size() / (double)blockSize);

				for (int idx = 0; idx < numBatches; idx++) {
					BatchData data;

					int begin = idx * blockSize;
					int end = begin + blockSize;
					if (end > pb->particles.size()) {
						end = pb->particles.size();
					}

					data.particlesBegin = pb->particles.begin() + begin;
					data.particlesEnd = pb->particles.begin() + end;

					data.refsBegin = pb->refs.begin() + begin;
					data.refsEnd = pb->refs.begin() + end;

					data.ufieldView = ugridView;
					data.vfieldView = vgridView;
					data.wfieldView = wgridView;

					data.batchOffset = bindex;
					data.indexOffset = indexOffset;
					data.positionOffset = Vector3f(positionOffset);

					batchData.push_back(data);
				}
			}
		}
	}
}

void ParticleAdvector::tricubicInterpolateBatch(std::vector<BatchData>& batches, std::vector<Vector3f>& output) {
	// prepare inputs
	std::vector<Vector3f> positionData;
	positionData.reserve(maxThreadsPerBlock * batches.size());

	for (int i = 0; i < batches.size(); i++) {
		BatchData b = batches[i];

		int numPoints = b.particlesEnd - b.particlesBegin;
		int numPadding = maxThreadsPerBlock - numPoints;
		Vector3f defaultPos(b.positionOffset);

		positionData.insert(positionData.end(), b.particlesBegin, b.particlesEnd);
		for (int j = 0; j < numPadding; j++) {
			positionData.push_back(defaultPos);
		}
	}

	std::vector<float> vfieldData;
	int vdataSize = 3 * sizeof(float) * (batchDim + 3) * (batchDim + 4) * (batchDim + 4);
	int numElems = vdataSize / sizeof(float);
	vfieldData.reserve(numElems);

	for (int idx = 0; idx < batches.size(); idx++) {
		BatchData b = batches[idx];
		for (int k = 0; k < b.ufieldView.depth; k++) {
			for (int j = 0; j < b.ufieldView.height; j++) {
				for (int i = 0; i < b.ufieldView.width; i++) {
					vfieldData.push_back(b.ufieldView(i, j, k));
				}
			}
		}

		for (int k = 0; k < b.vfieldView.depth; k++) {
			for (int j = 0; j < b.vfieldView.height; j++) {
				for (int i = 0; i < b.vfieldView.width; i++) {
					vfieldData.push_back(b.vfieldView(i, j, k));
				}
			}
		}

		for (int k = 0; k < b.wfieldView.depth; k++) {
			for (int j = 0; j < b.wfieldView.height; j++) {
				for (int i = 0; i < b.wfieldView.width; i++) {
					vfieldData.push_back(b.wfieldView(i, j, k));
				}
			}
		}
	}

	std::vector<Point3i> offsetData;
	offsetData.reserve(batches.size());
	for (int i = 0; i < batches.size(); i++) {
		offsetData.push_back(batches[i].batchOffset);
	}

	// launch kernel
	std::vector<Vector3f> batchOutput;
	cuda::launch_tricubicinterpolate_kernel(maxThreadsPerBlock, batches.size(), positionData, vfieldData, offsetData, cellsize, batchOutput);

	// copy data back into output
	for (int i = 0; i < batches.size(); i++) {
		int offset = i * maxThreadsPerBlock;
		int dataoffset = 0;

		BatchData batch = batches[i];
		auto begin = batch.refsBegin;
		auto end = batch.refsEnd;
		for (auto it = begin; it != end; it++) {
			output[*it] = batchOutput[offset + dataoffset];
			dataoffset++;
		}
	}
}

} // namespace foc