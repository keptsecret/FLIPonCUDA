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

} // namespace foc