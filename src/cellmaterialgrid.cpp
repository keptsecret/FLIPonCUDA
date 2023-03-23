#include "cellmaterialgrid.h"

namespace foc {

CellMaterialGrid::CellMaterialGrid() :
		isize(0), jsize(0), ksize(0) {
}

CellMaterialGrid::CellMaterialGrid(int i, int j, int k) :
		isize(i), jsize(j), ksize(k), grid(i, j, k, Material::air) {
}

Material CellMaterialGrid::operator()(int i, int j, int k) {
	return grid(i, j, k);
}

void CellMaterialGrid::fill(Material m) {
	grid.fill(m);
}

void CellMaterialGrid::set(int i, int j, int k, Material m) {
	grid.set(i, j, k, m);
}

bool CellMaterialGrid::isCellAir(int i, int j, int k) {
	return grid(i, j, k) == Material::air;
}

bool CellMaterialGrid::isCellFluid(int i, int j, int k) {
	return grid(i, j, k) == Material::fluid;
}

bool CellMaterialGrid::isCellSolid(int i, int j, int k) {
	return grid(i, j, k) == Material::solid;
}

bool CellMaterialGrid::isFaceBorderingMaterialU(int i, int j, int k, Material m) {
	if (i == grid.width) {
		return grid(i - 1, j, k) == m;
	} else if (i > 0) {
		return grid(i, j, k) == m || grid(i - 1, j, k) == m;
	} else {
		return grid(i, j, k) == m;
	}
}

bool CellMaterialGrid::isFaceBorderingMaterialV(int i, int j, int k, Material m) {
	if (j == grid.height) {
		return grid(i, j - 1, k) == m;
	} else if (j > 0) {
		return grid(i, j, k) == m || grid(i, j - 1, k) == m;
	} else {
		return grid(i, j, k) == m;
	}
}

bool CellMaterialGrid::isFaceBorderingMaterialW(int i, int j, int k, Material m) {
	if (k == grid.depth) {
		return grid(i, j, k - 1) == m;
	} else if (k > 0) {
		return grid(i, j, k) == m || grid(i, j, k - 1) == m;
	} else {
		return grid(i, j, k) == m;
	}
}

} // namespace foc