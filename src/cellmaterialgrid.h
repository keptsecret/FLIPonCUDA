#ifndef FOC_CELLMATERIALGRID_H
#define FOC_CELLMATERIALGRID_H

#include "foc.h"
#include "array3d.h"

namespace foc {

enum class Material : char {
	air = 0,
	fluid = 1,
	solid = 2
};

class CellMaterialGrid {
public:
	CellMaterialGrid();

	CellMaterialGrid(int i, int j, int k);

	Material operator()(int i, int j, int k);

	void fill(Material m);
	void set(int i, int j, int k, Material m);

	bool isCellAir(int i, int j, int k);
	bool isCellFluid(int i, int j, int k);
	bool isCellSolid(int i, int j, int k);

	int isize, jsize, ksize;

private:
	Array3D<Material> grid;
};

} // namespace foc

#endif // FOC_CELLMATERIALGRID_H
