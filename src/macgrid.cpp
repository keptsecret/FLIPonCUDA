#include "macgrid.h"

namespace foc {

MACGrid::MACGrid() :
		isize(10), jsize(10), ksize(10), cellsize(0.1) {
	initializeVelocityGrids();
}

MACGrid::MACGrid(int width, int height, int depth, double cellsize) :
		isize(width), jsize(height), ksize(depth), cellsize(cellsize) {
	initializeVelocityGrids();
}

void MACGrid::initializeVelocityGrids() {
	u = Array3D<float>(isize + 1, jsize, ksize, 0.0f);
	v = Array3D<float>(isize, jsize + 1, ksize, 0.0f);
	w = Array3D<float>(isize, jsize, ksize + 1, 0.0f);

	u.setOutOfRangeValue(0.0);
	v.setOutOfRangeValue(0.0);
	w.setOutOfRangeValue(0.0);
}

void MACGrid::clearU() {
	u.fill(0.0f);
}

void MACGrid::clearV() {
	v.fill(0.0f);
}

void MACGrid::clearW() {
	w.fill(0.0f);
}

void MACGrid::clear() {
	clearU();
	clearV();
	clearW();
}

Array3D<float>* MACGrid::getArrayU() {
	return &u;
}

Array3D<float>* MACGrid::getArrayV() {
	return &v;
}

Array3D<float>* MACGrid::getArrayW() {
	return &w;
}

float* MACGrid::getRawArrayU() {
	return u.data();
}

float* MACGrid::getRawArrayV() {
	return v.data();
}

float* MACGrid::getRawArrayW() {
	return w.data();
}

float MACGrid::U(int i, int j, int k) {
	if (!isGridIndexInRange(i, j, k, isize + 1, jsize, ksize)) {
		return OUT_OF_RANGE_VALUE;
	}

	return u(i, j, k);
}

float MACGrid::V(int i, int j, int k) {
	if (!isGridIndexInRange(i, j, k, isize, jsize + 1, ksize)) {
		return OUT_OF_RANGE_VALUE;
	}

	return v(i, j, k);
}

float MACGrid::W(int i, int j, int k) {
	if (!isGridIndexInRange(i, j, k, isize, jsize, ksize + 1)) {
		return OUT_OF_RANGE_VALUE;
	}

	return w(i, j, k);
}

void MACGrid::setU(int i, int j, int k, double val) {
	u.set(i, j, k, val);
}

void MACGrid::setV(int i, int j, int k, double val) {
	v.set(i, j, k, val);
}

void MACGrid::setW(int i, int j, int k, double val) {
	w.set(i, j, k, val);
}

Vector3f MACGrid::getVelocityCell(int i, int j, int k) {
	double xa = 0.5 * (U(i + 1, j, k) + U(i, j, k));
	double ya = 0.5 * (V(i, j + 1, k) + V(i, j, k));
	double za = 0.5 * (W(i, j, k + 1) + W(i, j, k));

	return Vector3f(xa, ya, za);
}

float MACGrid::getVelocityMagCell(int i, int j, int k) {
	double magSq = getVelocityMagSqCell(i, j, k);
	if (magSq > 0.0) {
		return static_cast<float>(sqrt(magSq));
	}

	return 0.f;
}

float MACGrid::getVelocityMagSqCell(int i, int j, int k) {
	double xa = 0.5 * (U(i + 1, j, k) + U(i, j, k));
	double ya = 0.5 * (V(i, j + 1, k) + V(i, j, k));
	double za = 0.5 * (W(i, j, k + 1) + W(i, j, k));

	return static_cast<float>(xa * xa + ya * ya + za * za);
}

float MACGrid::getMaxVelocityMag() {
	double maxSq = 0.0;
	for (int k = 0; k < ksize; k++) {
		for (int j = 0; j < jsize; j++) {
			for (int i = 0; i < isize; i++) {
				double m = getVelocityMagSqCell(i, j, k);
				maxSq = fmax(maxSq, m);
			}
		}
	}

	double max = maxSq;
	if (maxSq > 0.0) {
		max = sqrt(maxSq);
	}

	return static_cast<float>(max);
}

Vector3f MACGrid::getVelocityFaceU(int i, int j, int k) {
	// Shift reference coordinate to the left. The formula used is for calculating
	// u(i+1/2, j, k). If we keep original (i,j,k) coordinate, then using the formula
	// would calculate u(i+3/2, j, k) instead. The same will be done for the V and W
	// faces, shifting back in the respective direction.

	i--;

	double vx = U(i + 1, j, k);
	double vy = 0.25 * (V(i, j, k) + V(i, j + 1, k) + V(i + 1, j, k) + V(i + 1, j + 1, k));
	double vz = 0.25 * (W(i, j, k) + W(i, j, k + 1) + W(i + 1, j, k) + W(i + 1, j, k + 1));

	return Vector3f(vx, vy, vz);
}

Vector3f MACGrid::getVelocityFaceV(int i, int j, int k) {
	j--;

	double vx = 0.25 * (U(i, j, k) + U(i + 1, j, k) + U(i, j + 1, k) + U(i + 1, j + 1, k));
	double vy = V(i, j + 1, k);
	double vz = 0.25 * (W(i, j, k) + W(i, j, k + 1) + W(i, j + 1, k) + W(i, j + 1, k + 1));

	return Vector3f(vx, vy, vz);
}

Vector3f MACGrid::getVelocityFaceW(int i, int j, int k) {
	k--;

	double vx = 0.25 * (U(i, j, k) + U(i + 1, j, k) + U(i, j, k + 1) + U(i + 1, j, k + 1));
	double vy = 0.25 * (V(i, j, k) + V(i, j + 1, k) + V(i, j, k + 1) + V(i, j + 1, k + 1));
	double vz = W(i, j, k + 1);

	return Vector3f(vx, vy, vz);
}

double MACGrid::interpolateU(double x, double y, double z) {
	if (!isPositionInGrid(x, y, z, cellsize, isize, jsize, ksize)) {
		return 0.0;
	}

	y -= 0.5 * cellsize;
	z -= 0.5 * cellsize;

	Point3i idx = positionToGridIndex(x, y, z, cellsize);
	Point3f pos = gridIndexToPosition(idx.x, idx.y, idx.z, cellsize);

	double invcs = 1.0 / cellsize;
	double ix = (x - pos.x) * invcs;
	double iy = (y - pos.y) * invcs;
	double iz = (z - pos.z) * invcs;

	int refi = idx.x - 1;
	int refj = idx.y - 1;
	int refk = idx.z - 1;

	double points[4][4][4];
	for (int pk = 0; pk < 4; pk++) {
		for (int pj = 0; pj < 4; pj++) {
			for (int pi = 0; pi < 4; pi++) {
				points[pk][pj][pi] = U(pi + refi, pj + refj, pk + refk);
			}
		}
	}

	return tricubicInterpolate(points, ix, iy, iz);
}

double MACGrid::interpolateV(double x, double y, double z) {
	if (!isPositionInGrid(x, y, z, cellsize, isize, jsize, ksize)) {
		return 0.0;
	}

	x -= 0.5 * cellsize;
	z -= 0.5 * cellsize;

	Point3i idx = positionToGridIndex(x, y, z, cellsize);
	Point3f pos = gridIndexToPosition(idx.x, idx.y, idx.z, cellsize);

	double invcs = 1.0 / cellsize;
	double ix = (x - pos.x) * invcs;
	double iy = (y - pos.y) * invcs;
	double iz = (z - pos.z) * invcs;

	int refi = idx.x - 1;
	int refj = idx.y - 1;
	int refk = idx.z - 1;

	double points[4][4][4];
	for (int pk = 0; pk < 4; pk++) {
		for (int pj = 0; pj < 4; pj++) {
			for (int pi = 0; pi < 4; pi++) {
				points[pk][pj][pi] = V(pi + refi, pj + refj, pk + refk);
			}
		}
	}

	return tricubicInterpolate(points, ix, iy, iz);
}

double MACGrid::interpolateW(double x, double y, double z) {
	if (!isPositionInGrid(x, y, z, cellsize, isize, jsize, ksize)) {
		return 0.0;
	}

	x -= 0.5 * cellsize;
	y -= 0.5 * cellsize;

	Point3i idx = positionToGridIndex(x, y, z, cellsize);
	Point3f pos = gridIndexToPosition(idx.x, idx.y, idx.z, cellsize);

	double invcs = 1.0 / cellsize;
	double ix = (x - pos.x) * invcs;
	double iy = (y - pos.y) * invcs;
	double iz = (z - pos.z) * invcs;

	int refi = idx.x - 1;
	int refj = idx.y - 1;
	int refk = idx.z - 1;

	double points[4][4][4];
	for (int pk = 0; pk < 4; pk++) {
		for (int pj = 0; pj < 4; pj++) {
			for (int pi = 0; pi < 4; pi++) {
				points[pk][pj][pi] = W(pi + refi, pj + refj, pk + refk);
			}
		}
	}

	return tricubicInterpolate(points, ix, iy, iz);
}

Vector3f MACGrid::getVelocityAt(double x, double y, double z) {
	double xvel = interpolateU(x, y, z);
	double yvel = interpolateV(x, y, z);
	double zvel = interpolateW(x, y, z);

	return Vector3f(xvel, yvel, zvel);
}

Vector3f MACGrid::getVelocityAt(Point3f pos) {
	return getVelocityAt(pos.x, pos.y, pos.z);
}

void MACGrid::getGridDimensions(int* i, int* j, int* k) {
	*i = isize;
	*j = jsize;
	*k = ksize;
}

void MACGrid::resetExtrapolatedFluidVelocities(CellMaterialGrid& materialGrid) {
	for (int k = 0; k < ksize; k++) {
		for (int j = 0; j < jsize; j++) {
			for (int i = 0; i < isize + 1; i++) {
				if (!materialGrid.isFaceBorderingMaterialU(i, j, k, Material::fluid)) {
					setU(i, j, k, 0.0);
				}
			}
		}
	}

	for (int k = 0; k < ksize; k++) {
		for (int j = 0; j < jsize + 1; j++) {
			for (int i = 0; i < isize; i++) {
				if (!materialGrid.isFaceBorderingMaterialV(i, j, k, Material::fluid)) {
					setV(i, j, k, 0.0);
				}
			}
		}
	}

	for (int k = 0; k < ksize + 1; k++) {
		for (int j = 0; j < jsize; j++) {
			for (int i = 0; i < isize; i++) {
				if (!materialGrid.isFaceBorderingMaterialW(i, j, k, Material::fluid)) {
					setW(i, j, k, 0.0);
				}
			}
		}
	}
}

void MACGrid::updateExtrapolationLayers(CellMaterialGrid& materialGrid, Array3D<int>& layerGrid) {
	for (int k = 0; k < materialGrid.ksize; k++) {
		for (int j = 0; j < materialGrid.jsize; j++) {
			for (int i = 0; i < materialGrid.isize; i++) {
				if (materialGrid.isCellFluid(i, j, k)) {
					layerGrid.set(i, j, k, 0);
				}
			}
		}
	}

	for (int layer = 1; layer <= numExtrapolationLayers; layer++) {
		Point3i neighbors[6];

		for (int k = 0; k < layerGrid.depth; k++) {
			for (int j = 0; j < layerGrid.height; j++) {
				for (int i = 0; i < layerGrid.width; i++) {
					if (layerGrid(i, j, k) == layer - 1 && materialGrid.isCellSolid(i, j, k)) {
						getNeighborGridIndices6(i, j, k, neighbors);
						for (int idx = 0; idx < 6; idx++) {
							Point3i nb = neighbors[i];

							if (isGridIndexInRange(nb, isize, jsize, ksize) && layerGrid(nb.x, nb.y, nb.z) == -1 && materialGrid.isCellSolid(nb.x, nb.y, nb.z)) {
								layerGrid.set(nb.x, nb.y, nb.z, layer);
							}
						}
					}
				}
			}
		}
	}
}

double MACGrid::getExtrapolatedVelocityForFaceU(int i, int j, int k, int layerIdx, Array3D<int>& layerGrid) {
	Point3i neighbors[6];
	getNeighborGridIndices6(i, j, k, neighbors);

	double sum = 0.0, weights = 0.0;
	for (int idx = 0; idx < 6; idx++) {
		Point3i nb = neighbors[idx];
		if (isGridIndexInRange(i, j, k, isize + 1, jsize, ksize) && isFaceBorderingGridValueU(nb.x, nb.y, nb.z, layerIdx - 1, layerGrid)) {
			sum += U(nb.x, nb.y, nb.z);
			weights++;
		}
	}

	if (sum == 0.0) {
		return 0.0;
	}

	return sum / weights;
}

double MACGrid::getExtrapolatedVelocityForFaceV(int i, int j, int k, int layerIdx, Array3D<int>& layerGrid) {
	Point3i neighbors[6];
	getNeighborGridIndices6(i, j, k, neighbors);

	double sum = 0.0, weights = 0.0;
	for (int idx = 0; idx < 6; idx++) {
		Point3i nb = neighbors[idx];
		if (isGridIndexInRange(i, j, k, isize, jsize + 1, ksize) && isFaceBorderingGridValueV(nb.x, nb.y, nb.z, layerIdx - 1, layerGrid)) {
			sum += V(nb.x, nb.y, nb.z);
			weights++;
		}
	}

	if (sum == 0.0) {
		return 0.0;
	}

	return sum / weights;
}

double MACGrid::getExtrapolatedVelocityForFaceW(int i, int j, int k, int layerIdx, Array3D<int>& layerGrid) {
	Point3i neighbors[6];
	getNeighborGridIndices6(i, j, k, neighbors);

	double sum = 0.0, weights = 0.0;
	for (int idx = 0; idx < 6; idx++) {
		Point3i nb = neighbors[idx];
		if (isGridIndexInRange(i, j, k, isize, jsize, ksize + 1) && isFaceBorderingGridValueW(nb.x, nb.y, nb.z, layerIdx - 1, layerGrid)) {
			sum += W(nb.x, nb.y, nb.z);
			weights++;
		}
	}

	if (sum == 0.0) {
		return 0.0;
	}

	return sum / weights;
}

void MACGrid::extrapolateVelocitiesForLayerIndexU(int idx, CellMaterialGrid& materialGrid, Array3D<int>& layerGrid) {
	for (int k = 0; k < ksize; k++) {
		for (int j = 0; j < jsize; j++) {
			for (int i = 0; i < isize + 1; i++) {
				bool isExtrapolated = isFaceBorderingGridValueU(i, j, k, idx, layerGrid) && isFaceBorderingGridValueU(i - 1, j, k, idx, layerGrid) && (!materialGrid.isFaceBorderingMaterialU(i, j, k, Material::solid));
				if (isExtrapolated) {
					double vel = getExtrapolatedVelocityForFaceU(i, j, k, idx, layerGrid);
					setU(i, j, k, vel);
				}
			}
		}
	}
}

void MACGrid::extrapolateVelocitiesForLayerIndexV(int idx, CellMaterialGrid& materialGrid, Array3D<int>& layerGrid) {
	for (int k = 0; k < ksize; k++) {
		for (int j = 0; j < jsize + 1; j++) {
			for (int i = 0; i < isize; i++) {
				bool isExtrapolated = isFaceBorderingGridValueV(i, j, k, idx, layerGrid) && isFaceBorderingGridValueV(i, j - 1, k, idx, layerGrid) && (!materialGrid.isFaceBorderingMaterialV(i, j, k, Material::solid));
				if (isExtrapolated) {
					double vel = getExtrapolatedVelocityForFaceV(i, j, k, idx, layerGrid);
					setV(i, j, k, vel);
				}
			}
		}
	}
}

void MACGrid::extrapolateVelocitiesForLayerIndexW(int idx, CellMaterialGrid& materialGrid, Array3D<int>& layerGrid) {
	for (int k = 0; k < ksize + 1; k++) {
		for (int j = 0; j < jsize; j++) {
			for (int i = 0; i < isize; i++) {
				bool isExtrapolated = isFaceBorderingGridValueW(i, j, k, idx, layerGrid) && isFaceBorderingGridValueW(i, j, k - 1, idx, layerGrid) && (!materialGrid.isFaceBorderingMaterialW(i, j, k, Material::solid));
				if (isExtrapolated) {
					double vel = getExtrapolatedVelocityForFaceW(i, j, k, idx, layerGrid);
					setW(i, j, k, vel);
				}
			}
		}
	}
}

void MACGrid::extrapolateVelocitiesForLayerIndex(int idx, CellMaterialGrid& materialGrid, Array3D<int>& layerGrid) {
	extrapolateVelocitiesForLayerIndexU(idx, materialGrid, layerGrid);
	extrapolateVelocitiesForLayerIndexV(idx, materialGrid, layerGrid);
	extrapolateVelocitiesForLayerIndexW(idx, materialGrid, layerGrid);
}

void MACGrid::extrapolateVelocityField(CellMaterialGrid& materialGrid, int numLayers) {
	numExtrapolationLayers = numLayers;

	Array3D<int> layerGrid(isize, jsize, ksize, -1);

	resetExtrapolatedFluidVelocities(materialGrid);
	updateExtrapolationLayers(materialGrid, layerGrid);

	for (int i = 1; i <= numLayers; i++) {
		extrapolateVelocitiesForLayerIndex(i, materialGrid, layerGrid);
	}
}

} // namespace foc