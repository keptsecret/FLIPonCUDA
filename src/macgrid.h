#ifndef FOC_MACGRID_H
#define FOC_MACGRID_H

#include "array3d.h"
#include "cellmaterialgrid.h"
#include "foc.h"

namespace foc {

class MACGrid {
public:
	MACGrid();

	MACGrid(int width, int height, int depth, double cellsize);

	~MACGrid() {}

	void clearU();
	void clearV();
	void clearW();
	void clear();

	Array3D<float>* getArrayU();
	Array3D<float>* getArrayV();
	Array3D<float>* getArrayW();

	float* getRawArrayU();
	float* getRawArrayV();
	float* getRawArrayW();

	float U(int i, int j, int k);
	float V(int i, int j, int k);
	float W(int i, int j, int k);

	void setU(int i, int j, int k, double val);
	void setV(int i, int j, int k, double val);
	void setW(int i, int j, int k, double val);

	Vector3f getVelocityCell(int i, int j, int k);
	float getVelocityMagCell(int i, int j, int k);
	float getVelocityMagSqCell(int i, int j, int k);
	float getMaxVelocityMag();

	Vector3f getVelocityFaceU(int i, int j, int k);
	Vector3f getVelocityFaceV(int i, int j, int k);
	Vector3f getVelocityFaceW(int i, int j, int k);

	Vector3f getVelocityAt(double x, double y, double z);
	Vector3f getVelocityAt(Point3f pos);

	void extrapolateVelocityField(CellMaterialGrid& materialGrid, int numLayers);

private:
	void initializeVelocityGrids();

	void resetExtrapolatedFluidVelocities(CellMaterialGrid& materialGrid);
	void updateExtrapolationLayers(CellMaterialGrid& materialGrid, Array3D<int>& layerGrid);

	double getExtrapolatedVelocityForFaceU(int i, int j, int k, int layerIdx, Array3D<int>& layerGrid);
	double getExtrapolatedVelocityForFaceV(int i, int j, int k, int layerIdx, Array3D<int>& layerGrid);
	double getExtrapolatedVelocityForFaceW(int i, int j, int k, int layerIdx, Array3D<int>& layerGrid);

	void extrapolateVelocitiesForLayerIndexU(int idx, CellMaterialGrid& materialGrid, Array3D<int>& layerGrid);
	void extrapolateVelocitiesForLayerIndexV(int idx, CellMaterialGrid& materialGrid, Array3D<int>& layerGrid);
	void extrapolateVelocitiesForLayerIndexW(int idx, CellMaterialGrid& materialGrid, Array3D<int>& layerGrid);
	void extrapolateVelocitiesForLayerIndex(int idx, CellMaterialGrid& materialGrid, Array3D<int>& layerGrid);

	template <typename T>
	bool isFaceBorderingGridValueU(int i, int j, int k, int value, Array3D<T>& grid) {
		if (i == grid.width) {
			return grid(i - 1, j, k) == value;
		} else if (i > 0) {
			return grid(i - 1, j, k) == value || grid(i, j, k) == value;
		} else {
			return grid(i, j, k) == value;
		}
	}

	template <typename T>
	bool isFaceBorderingGridValueV(int i, int j, int k, int value, Array3D<T>& grid) {
		if (j == grid.height) {
			return grid(i, j - 1, k) == value;
		} else if (j > 0) {
			return grid(i, j - 1, k) == value || grid(i, j, k) == value;
		} else {
			return grid(i, j, k) == value;
		}
	}

	template <typename T>
	bool isFaceBorderingGridValueW(int i, int j, int k, int value, Array3D<T>& grid) {
		if (k == grid.depth) {
			return grid(i, j, k - 1) == value;
		} else if (k > 0) {
			return grid(i, j, k - 1) == value || grid(i, j, k) == value;
		} else {
			return grid(i, j, k) == value;
		}
	}

private:
	int numExtrapolationLayers;

	int isize, jsize, ksize;
	double dcell;
	Array3D<float> u, v, w;
};

} // namespace foc

#endif // FOC_MACGRID_H
