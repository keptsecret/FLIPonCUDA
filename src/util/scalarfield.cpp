#include "scalarfield.h"

#ifdef FOC_BUILD_GPU
#include <gpu/kernels.h>
#endif

#include <algorithm>

namespace foc {

ScalarField::ScalarField() {
}

ScalarField::ScalarField(int i, int j, int k, double dx) :
		isize(i), jsize(j), ksize(k), cellsize(dx), field(i, j, k, false), isVertexSolid(i, j, k, false), isVertexSet(i, j, k, false) {
}

void ScalarField::setPointRadius(double r) {
	radius = r;
	invRadius = 1 / r;
	coef1 = (4.0 / 9.0) * (1.0 / (r * r * r * r * r * r));
	coef2 = (17.0 / 9.0) * (1.0 / (r * r * r * r));
	coef3 = (22.0 / 9.0) * (1.0 / (r * r));
}

void ScalarField::setFieldOffset(Vector3f offset) {
	gridOffset = offset;
}

void ScalarField::enableWeightField() {
	if (isWeightFieldEnabled) {
		return;
	}

	weightField = Array3D<float>(isize, jsize, ksize, 0.0f);
	isWeightFieldEnabled = true;
}

void ScalarField::applyWeightField() {
	if (!isWeightFieldEnabled) {
		return;
	}

	for (int k = 0; k < ksize; k++) {
		for (int j = 0; j < jsize; j++) {
			for (int i = 0; i < isize; i++) {
				float weight = weightField(i, j, k);
				if (weight > 0.0) {
					float v = field(i, j, k) / weight;
					setScalarFieldValue(i, j, k, v);
				}
			}
		}
	}
}

void ScalarField::addPoint(Point3f p, double r) {
	setPointRadius(r);
	addPoint(p);
}

void ScalarField::addPoint(Point3f p) {
	p -= gridOffset;

	Point3i gmin, gmax;
	getGridIndexBounds(p, radius, cellsize, isize, jsize, ksize, &gmin, &gmax);

	Point3f gpos;
	Vector3f v;
	double rsq = radius * radius;
	double distsq;
	double weight;
	for (int k = gmin.z; k <= gmax.z; k++) {
		for (int j = gmin.y; j <= gmax.y; j++) {
			for (int i = gmin.x; i <= gmax.x; i++) {
				if (isMaxScalarFieldThresholdSet && field(i, j, k) > maxScalarFieldThreshold) {
					continue;
				}

				gpos = gridIndexToPosition(i, j, k, cellsize);
				v = gpos - p;
				distsq = dot(v, v);
				if (distsq < rsq) {
					weight = evaluateTricubicFieldFunctionForRadiusSquared(distsq);
					addScalarFieldValue(i, j, k, weight);

					if (isWeightFieldEnabled) {
						weightField.set(i, j, k, weightField(i, j, k) + (float)weight);
					}
				}
			}
		}
	}
}

void ScalarField::addCuboid(Point3f pmin, Point3f pmax) {
	pmin -= gridOffset;

	Point3i gmin, gmax;
	Bounds3f bbox = Bounds3f(pmin, pmax);
	getGridIndexBounds(bbox, cellsize, isize, jsize, ksize, &gmin, &gmax);

	Point3f gpos;
	double eps = 10e-6;
	for (int k = gmin.z; k <= gmax.z; k++) {
		for (int j = gmin.y; j <= gmax.y; j++) {
			for (int i = gmin.x; i <= gmax.x; i++) {
				if (isMaxScalarFieldThresholdSet && field(i, j, k) > maxScalarFieldThreshold) {
					continue;
				}

				gpos = gridIndexToPosition(i, j, k, cellsize);
				if (bInside(gpos, bbox)) {
					addScalarFieldValue(i, j, k, surfaceThreshold + eps);

					if (isWeightFieldEnabled) {
						weightField.set(i, j, k, (float)(weightField(i, j, k) + surfaceThreshold + eps));
					}
				}
			}
		}
	}
}

void ScalarField::addPointValue(Point3f p, double val) {
	p -= gridOffset;

	Point3i gmin, gmax;
	getGridIndexBounds(p, radius, cellsize, isize, jsize, ksize, &gmin, &gmax);

	Point3f gpos;
	Vector3f v;
	double rsq = radius * radius;
	double distsq;
	double weight;
	for (int k = gmin.z; k <= gmax.z; k++) {
		for (int j = gmin.y; j <= gmax.y; j++) {
			for (int i = gmin.x; i <= gmax.x; i++) {
				if (isMaxScalarFieldThresholdSet && field(i, j, k) > maxScalarFieldThreshold) {
					continue;
				}

				gpos = gridIndexToPosition(i, j, k, cellsize);
				v = gpos - p;
				distsq = dot(v, v);
				if (distsq < rsq) {
					weight = evaluateTricubicFieldFunctionForRadiusSquared(distsq);
					addScalarFieldValue(i, j, k, weight * val);

					if (isWeightFieldEnabled) {
						weightField.set(i, j, k, weightField(i, j, k) + (float)weight);
					}
				}
			}
		}
	}
}

void ScalarField::addPoints(std::vector<Vector3f>& points) {
	if (!field.isOutOfRangeValueSet()) {
		field.setOutOfRangeValue(0.0);
	}

	std::vector<PointData> pointData;
	initializePointData(points, pointData);

	Point3i batchDataDims(std::ceil((double)isize / (double)batchDim), std::ceil((double)jsize / (double)batchDim), std::ceil((double)ksize / (double)batchDim));
	Array3D<BatchData> batchDataGrid(batchDataDims.x, batchDataDims.y, batchDataDims.z);
	initializeBatchDataGrid(pointData, batchDataGrid);

	std::vector<BatchDataRef> batchDataQueue;
	initializeBatchRefs(batchDataGrid, batchDataQueue);

	int batchSize = getMaxBatchesPerPointDataComputation();

	while (!batchDataQueue.empty()) {
		updateBatchMinimumValues(batchDataGrid);

		std::vector<BatchDataRef> batches;
		getNextBatches(batchDataQueue, batchDataGrid, batches, batchSize);
		computePointScalarField(batches, batchDataGrid);
	}

	if (!field.isOutOfRangeValueSet()) {
		field.setOutOfRangeValue();
	}
}

void ScalarField::addPointValues(std::vector<Vector3f>& points, std::vector<float>& values) {
	addPointValues(points, values, radius, gridOffset, cellsize);
}

void ScalarField::addPointValues(std::vector<Vector3f>& points, std::vector<float>& values, double radius, Vector3f offset, double dx) {
	if (!field.isOutOfRangeValueSet()) {
		field.setOutOfRangeValue(0.0);
	}
	if (isWeightFieldEnabled && !weightField.isOutOfRangeValueSet()) {
		weightField.setOutOfRangeValue(0.0);
	}

	std::vector<PointData> pointData;
	initializePointData(points, values, pointData);

	Point3i batchDataDims(std::ceil((double)isize / (double)batchDim), std::ceil((double)jsize / (double)batchDim), std::ceil((double)ksize / (double)batchDim));
	Array3D<BatchData> batchDataGrid(batchDataDims.x, batchDataDims.y, batchDataDims.z);
	initializeBatchDataGrid(pointData, batchDataGrid);

	std::vector<BatchDataRef> batchDataQueue;
	initializeBatchRefs(batchDataGrid, batchDataQueue);

	int batchSize = getMaxBatchesPerPointDataComputation();

	while (!batchDataQueue.empty()) {
		updateBatchMinimumValues(batchDataGrid);

		std::vector<BatchDataRef> batches;
		getNextBatches(batchDataQueue, batchDataGrid, batches, batchSize);
		if (!isWeightFieldEnabled) {
			computePointValueScalarField(batches, batchDataGrid);
		} else {
			computePointValueScalarWeightField(batches, batchDataGrid);
		}
	}

	if (!field.isOutOfRangeValueSet()) {
		field.setOutOfRangeValue();
	}
	if (isWeightFieldEnabled && !weightField.isOutOfRangeValueSet()) {
		weightField.setOutOfRangeValue();
	}
}

void ScalarField::setScalarFieldValue(int i, int j, int k, double value) {
	field.set(i, j, k, value);
	isVertexSet.set(i, j, k, true);
}

void ScalarField::addScalarFieldValue(int i, int j, int k, double value) {
	field.set(i, j, k, field(i, j, k) + value);
	isVertexSet.set(i, j, k, true);
}

double ScalarField::getScalarFieldValue(int i, int j, int k) {
	double val = field(i, j, k);
	if (isVertexSolid(i, j, k) && val > surfaceThreshold) {
		val = surfaceThreshold;
	}

	return val;
}

Array3D<float>* ScalarField::getScalarFieldPointer() {
	return &field;
}

Array3D<float>* ScalarField::getWeightFieldPointer() {
	return &weightField;
}

void ScalarField::setMaterialGrid(CellMaterialGrid& materialGrid) {
	Point3i vertices[8];
	for (int k = 0; k < ksize - 1; k++) {
		for (int j = 0; j < jsize - 1; j++) {
			for (int i = 0; i < isize - 1; i++) {
				if (materialGrid.isCellSolid(i, j, k)) {
					getGridIndexVertices(i, j, k, vertices);
					for (int idx = 0; idx < 8; idx++) {
						isVertexSolid.set(vertices[idx].x, vertices[idx].y, vertices[idx].z, true);
					}
				}
			}
		}
	}
}

double ScalarField::getScalarFieldValueAtCellCenter(int i, int j, int k) {
	double sum = 0.0;
	sum += getScalarFieldValue(i, j, k);
	sum += getScalarFieldValue(i + 1, j, k);
	sum += getScalarFieldValue(i, j + 1, k);
	sum += getScalarFieldValue(i + 1, j + 1, k);
	sum += getScalarFieldValue(i, j, k + 1);
	sum += getScalarFieldValue(i + 1, j, k + 1);
	sum += getScalarFieldValue(i, j + 1, k + 1);
	sum += getScalarFieldValue(i + 1, j + 1, k + 1);

	return 0.125 * sum;
}

double ScalarField::evaluateTricubicFieldFunctionForRadiusSquared(double rsq) {
	return 1.0 - coef1 * rsq * rsq * rsq + coef2 * rsq * rsq - coef3 * rsq;
}

void ScalarField::initializePointData(std::vector<Vector3f>& points, std::vector<PointData>& pd) {
	float defaultValue = 0.0f;
	Vector3f offset = gridOffset - Vector3f(0.5 * cellsize, 0.5 * cellsize, 0.5 * cellsize);
	pd.reserve(points.size());

	for (int i = 0; i < points.size(); i++) {
		pd.push_back(PointData(points[i] - offset, defaultValue));
	}
}

void ScalarField::initializePointData(std::vector<Vector3f>& points, std::vector<float>& values, std::vector<PointData>& pd) {
	Vector3f offset = gridOffset - Vector3f(0.5 * cellsize, 0.5 * cellsize, 0.5 * cellsize);
	pd.reserve(points.size());

	for (int i = 0; i < points.size(); i++) {
		pd.push_back(PointData(points[i] - offset, values[i]));
	}
}

void ScalarField::initializeBatchDataGrid(std::vector<PointData>& points, Array3D<BatchData>& grid) {
	for (int k = 0; k < grid.depth; k++) {
		for (int j = 0; j < grid.height; j++) {
			for (int i = 0; i < grid.width; i++) {
				BatchData* batch = grid.getPointer(i, j, k);
				Point3i batchOffset(i, j, k);
				Point3i indexOffset(i * batchDim, j * batchDim, k * batchDim);

				Point3f positionOffset = gridIndexToPosition(indexOffset.x, indexOffset.y, indexOffset.z, cellsize);

				batch->fieldView = Array3DView<float>(batchDim, batchDim, batchDim, indexOffset.x, indexOffset.y, indexOffset.z, &field);
				if (isWeightFieldEnabled) {
					batch->weightfieldView = Array3DView<float>(batchDim, batchDim, batchDim, indexOffset.x, indexOffset.y, indexOffset.z, &weightField);
				}

				batch->batchOffset = batchOffset;
				batch->indexOffset = indexOffset;
				batch->positionOffset = Vector3f(positionOffset);
			}
		}
	}

	double batchdx = batchDim * cellsize;
	double batchdy = batchDim * cellsize;
	double batchdz = batchDim * cellsize;
	double invbatchdx = 1.0 / batchdx;
	double invbatchdy = 1.0 / batchdy;
	double invbatchdz = 1.0 / batchdz;

	Point3i gmax(grid.width, grid.height, grid.depth);

	for (int i = 0; i < points.size(); i++) {
		PointData pd = points[i];
		Vector3f p = pd.position;

		int bi = std::floor(p.x * invbatchdx);
		int bj = std::floor(p.y * invbatchdy);
		int bk = std::floor(p.z * invbatchdz);
		double bx = (double)bi * batchdx;
		double by = (double)bj * batchdy;
		double bz = (double)bk * batchdz;

		Bounds3f cbbox(Point3f(bx + radius, by + radius, bz + radius), Point3f(bx + batchdx - radius, by + batchdy - radius, bz + batchdz - radius));
		if (bInside(Point3f(p.x, p.y, p.z), cbbox) && isGridIndexInRange(bi, bj, bk, gmax.x, gmax.y, gmax.z)) {
			// sphere inside one grid cell
			BatchData* b = grid.getPointer(bi, bj, bk);
			b->particles.push_back(pd);
			continue;
		}

		// sphere inside at least two grid cells
		Vector3f minp(p.x - radius, p.y - radius, p.z - radius);
		Vector3f maxp(p.x + radius, p.y + radius, p.z + radius);
		int mini = fmax(std::floor(minp.x * invbatchdx), 0);
		int minj = fmax(std::floor(minp.y * invbatchdy), 0);
		int mink = fmax(std::floor(minp.z * invbatchdz), 0);
		int maxi = fmin(std::floor(maxp.x * invbatchdx), gmax.x - 1);
		int maxj = fmin(std::floor(maxp.y * invbatchdy), gmax.y - 1);
		int maxk = fmin(std::floor(maxp.z * invbatchdz), gmax.z - 1);

		for (int ck = mink; ck <= maxk; ck++) {
			for (int cj = minj; cj <= maxj; cj++) {
				for (int ci = mini; ci <= maxi; ci++) {
					BatchData* b = grid.getPointer(bi, bj, bk);
					b->particles.push_back(pd);
				}
			}
		}
	}
}

void ScalarField::initializeBatchRefs(Array3D<BatchData>& grid, std::vector<BatchDataRef>& queue) {
	for (int k = 0; k < grid.depth; k++) {
		for (int j = 0; j < grid.height; j++) {
			for (int i = 0; i < grid.width; i++) {
				BatchData* b = grid.getPointer(i, j, k);
				if (b->particles.size() == 0) {
					continue;
				}

				Point3i batchidx = b->batchOffset;
				int size = b->particles.size();
				int batchsize = maxParticlesPerBatch;
				for (int idx = 0; idx < size; idx += batchsize) {
					BatchDataRef ref;
					ref.batchDataIndex = batchidx;

					int begidx = idx;
					int endidx = begidx + batchsize;
					if (endidx > size) {
						endidx = size;
					}

					ref.particlesBegin = b->particles.begin() + begidx;
					ref.particlesEnd = b->particles.begin() + endidx;

					queue.push_back(ref);
				}
			}
		}
	}
	queue.shrink_to_fit();

	std::sort(queue.begin(), queue.end(), [](const BatchDataRef& b1, const BatchDataRef& b2) -> bool {
		return b1.particlesEnd - b1.particlesBegin < b2.particlesEnd - b2.particlesBegin;
	});
}

int ScalarField::getMaxBatchesPerPointDataComputation() {
	int pointDataSize = 4 * maxParticlesPerBatch * sizeof(float);
	int fieldDataSize = batchDim * batchDim * batchDim * sizeof(float);
	int offsetDataSize = 3 * sizeof(int);
	int totalSize = pointDataSize + fieldDataSize + offsetDataSize;

#ifndef FOC_BUILD_GPU
	fprintf(stderr, "Error:: ScalarField::getMaxBatchesPerPointDataComputation: this shouldn't be happening\n");
	return maxBatchesPerComputation;
#else
	size_t maxGlobalMem = cuda::getCUDADeviceGlobalMem();
	return fmin(std::floor((double)maxGlobalMem / (double)totalSize), maxBatchesPerComputation);
#endif
}

void ScalarField::updateBatchMinimumValues(Array3D<BatchData>& grid) {
	for (int k = 0; k < grid.depth; k++) {
		for (int j = 0; j < grid.height; j++) {
			for (int i = 0; i < grid.width; i++) {
				BatchData* b = grid.getPointer(i, j, k);
				float minval = MaxFloat;

				for (int vk = 0; vk < b->fieldView.depth; vk++) {
					for (int vj = 0; vj < b->fieldView.height; vj++) {
						for (int vi = 0; vi < b->fieldView.width; vi++) {
							if (b->fieldView.isIndexInParent(vi, vj, vk) && b->fieldView(vi, vj, vk) < minval) {
								minval = b->fieldView(vi, vj, vk);
							}
						}
					}
				}

				b->minFieldValue = minval;
			}
		}
	}
}

void ScalarField::getNextBatches(std::vector<BatchDataRef>& queue, Array3D<BatchData>& grid, std::vector<BatchDataRef>& batches, int n) {
	while (batches.size() < n && !queue.empty()) {
		BatchDataRef b = queue.back();
		queue.pop_back();

		if (isMaxScalarFieldThresholdSet) {
			BatchData* bd = grid.getPointer(b.batchDataIndex.x, b.batchDataIndex.y, b.batchDataIndex.z);
			float minval = bd->minFieldValue;

			if (minval < maxScalarFieldThreshold) {
				batches.push_back(b);
			}
		} else {
			batches.push_back(b);
		}
	}
}

void ScalarField::computePointScalarField(std::vector<BatchDataRef>& batches, Array3D<BatchData>& batchDataGrid) {
#ifdef FOC_BUILD_GPU
	int numParticles = getMaxNumParticlesInBatch(batches);

	std::vector<float> pointDataBuffer;
	fillPointDataBuffer(batches, batchDataGrid, numParticles, pointDataBuffer);

	std::vector<float> fieldDataBuffer;
	fillScalarFieldDataBuffer(batches, fieldDataBuffer);

	std::vector<Point3i> offsetDataBuffer;
	fillBatchOffsetDataBuffer(batches, offsetDataBuffer);

	cuda::launch_scalar_field_point_kernel(maxThreadsPerBlock, batches.size(), pointDataBuffer, fieldDataBuffer, offsetDataBuffer, numParticles, radius, cellsize);

	int bufferidx = 0;
	for (int bidx = 0; bidx < batches.size(); bidx++) {
		Point3i g = batches[bidx].batchDataIndex;
		Array3DView<float> fieldView = batchDataGrid(g.x, g.y, g.z).fieldView;

		for (int k = 0; k < fieldView.depth; k++) {
			for (int j = 0; j < fieldView.height; j++) {
				for (int i = 0; i < fieldView.width; i++) {
					fieldView.set(i, j, k, fieldView(i, j, k) + fieldDataBuffer[bufferidx]);
					bufferidx++;
				}
			}
		}
	}
#endif
}

void ScalarField::computePointValueScalarField(std::vector<BatchDataRef>& batches, Array3D<BatchData>& batchDataGrid) {
#ifdef FOC_BUILD_GPU
	int numParticles = getMaxNumParticlesInBatch(batches);

	std::vector<float> pointDataBuffer;
	fillPointValueDataBuffer(batches, batchDataGrid, numParticles, pointDataBuffer);

	std::vector<float> fieldDataBuffer;
	fillScalarFieldDataBuffer(batches, fieldDataBuffer);

	std::vector<Point3i> offsetDataBuffer;
	fillBatchOffsetDataBuffer(batches, offsetDataBuffer);

	cuda::launch_scalar_field_point_value_kernel(maxThreadsPerBlock, batches.size(), pointDataBuffer, fieldDataBuffer, offsetDataBuffer, numParticles, radius, cellsize);

	int bufferidx = 0;
	for (int bidx = 0; bidx < batches.size(); bidx++) {
		Point3i g = batches[bidx].batchDataIndex;
		Array3DView<float> fieldView = batchDataGrid(g.x, g.y, g.z).fieldView;

		for (int k = 0; k < fieldView.depth; k++) {
			for (int j = 0; j < fieldView.height; j++) {
				for (int i = 0; i < fieldView.width; i++) {
					fieldView.set(i, j, k, fieldView(i, j, k) + fieldDataBuffer[bufferidx]);
					bufferidx++;
				}
			}
		}
	}
#endif
}

void ScalarField::computePointValueScalarWeightField(std::vector<BatchDataRef>& batches, Array3D<BatchData>& batchDataGrid) {
#ifdef FOC_BUILD_GPU
	int numParticles = getMaxNumParticlesInBatch(batches);

	std::vector<float> pointDataBuffer;
	fillPointValueDataBuffer(batches, batchDataGrid, numParticles, pointDataBuffer);

	std::vector<float> fieldDataBuffer;
	fillScalarWeightFieldDataBuffer(batches, fieldDataBuffer);

	std::vector<Point3i> offsetDataBuffer;
	fillBatchOffsetDataBuffer(batches, offsetDataBuffer);

	cuda::launch_scalar_weight_field_point_value_kernel(maxThreadsPerBlock, batches.size(), pointDataBuffer, fieldDataBuffer, offsetDataBuffer, numParticles, radius, cellsize);

	int bufferidx = 0;
	int elementsPerBatch = batchDim * batchDim * batchDim;
	int weightfieldOffset = batches.size() * elementsPerBatch;
	for (int bidx = 0; bidx < batches.size(); bidx++) {
		Point3i g = batches[bidx].batchDataIndex;
		Array3DView<float> fieldView = batchDataGrid(g.x, g.y, g.z).fieldView;
		Array3DView<float> weightFieldView = batchDataGrid(g.x, g.y, g.z).weightfieldView;

		for (int k = 0; k < fieldView.depth; k++) {
			for (int j = 0; j < fieldView.height; j++) {
				for (int i = 0; i < fieldView.width; i++) {
					fieldView.set(i, j, k, fieldView(i, j, k) + fieldDataBuffer[bufferidx]);
					weightFieldView.set(i, j, k, weightFieldView(i, j, k) + fieldDataBuffer[bufferidx + weightfieldOffset]);
					bufferidx++;
				}
			}
		}
	}
#endif
}

int ScalarField::getMaxNumParticlesInBatch(std::vector<BatchDataRef>& batches) {
	int maxParticles = 0;
	for (int i = 0; i < batches.size(); i++) {
		int n = batches[i].particlesEnd - batches[i].particlesBegin;
		if (n > maxParticles) {
			maxParticles = n;
		}
	}

	return maxParticles;
}

void ScalarField::fillPointDataBuffer(std::vector<BatchDataRef>& batches, Array3D<BatchData>& grid, int numParticles, std::vector<float>& buffer) {
	int numElems = 3 * batches.size() * numParticles;
	buffer.reserve(numElems);

	Vector3f outOfRangePos(grid.width * batchDim * cellsize + 2 * radius,
			grid.height * batchDim * cellsize + 2 * radius,
			grid.depth * batchDim * cellsize + 2 * radius);

	for (int i = 0; i < batches.size(); i++) {
		BatchDataRef b = batches[i];

		int numPoints = b.particlesEnd - b.particlesBegin;
		int numPad = numParticles - numPoints;

		auto begin = b.particlesBegin;
		auto end = b.particlesEnd;
		for (auto it = begin; it != end; it++) {
			Vector3f p = (*it).position;
			buffer.push_back(p.x);
			buffer.push_back(p.y);
			buffer.push_back(p.z);
		}

		for (int j = 0; j < numPad; j++) {
			buffer.push_back(outOfRangePos.x);
			buffer.push_back(outOfRangePos.y);
			buffer.push_back(outOfRangePos.z);
		}
	}
}

void ScalarField::fillPointValueDataBuffer(std::vector<BatchDataRef>& batches, Array3D<BatchData>& grid, int numParticles, std::vector<float>& buffer) {
	int numElems = 4 * batches.size() * numParticles;
	buffer.reserve(numElems);

	Vector3f outOfRangePos(grid.width * batchDim * cellsize + 2 * radius,
			grid.height * batchDim * cellsize + 2 * radius,
			grid.depth * batchDim * cellsize + 2 * radius);

	float outOfRangeValue = 0.0f;

	for (int i = 0; i < batches.size(); i++) {
		BatchDataRef b = batches[i];

		int numPoints = b.particlesEnd - b.particlesBegin;
		int numPad = numParticles - numPoints;

		auto begin = b.particlesBegin;
		auto end = b.particlesEnd;
		for (auto it = begin; it != end; it++) {
			Vector3f p = (*it).position;
			float v = (*it).value;
			buffer.push_back(p.x);
			buffer.push_back(p.y);
			buffer.push_back(p.z);
			buffer.push_back(v);
		}

		for (int j = 0; j < numPad; j++) {
			buffer.push_back(outOfRangePos.x);
			buffer.push_back(outOfRangePos.y);
			buffer.push_back(outOfRangePos.z);
			buffer.push_back(outOfRangeValue);
		}
	}
}

void ScalarField::fillScalarFieldDataBuffer(std::vector<BatchDataRef>& batches, std::vector<float>& buffer) {
	int numElems = batches.size() * batchDim * batchDim * batchDim;
	buffer.reserve(numElems);
	for (int i = 0; i < numElems; i++) {
		buffer.push_back(0.0);
	}
}

void ScalarField::fillScalarWeightFieldDataBuffer(std::vector<BatchDataRef>& batches, std::vector<float>& buffer) {
	int numElems = 2 * batches.size() * batchDim * batchDim * batchDim;
	buffer.reserve(numElems);
	for (int i = 0; i < numElems; i++) {
		buffer.push_back(0.0);
	}
}

void ScalarField::fillBatchOffsetDataBuffer(std::vector<BatchDataRef>& batches, std::vector<Point3i>& buffer) {
	buffer.reserve(batches.size());
	for (int i = 0; i < batches.size(); i++) {
		buffer.push_back(batches[i].batchDataIndex);
	}
}

} // namespace foc