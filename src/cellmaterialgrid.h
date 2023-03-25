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

	bool isFaceBorderingMaterialU(int i, int j, int k, Material m);
	bool isFaceBorderingMaterialV(int i, int j, int k, Material m);
	bool isFaceBorderingMaterialW(int i, int j, int k, Material m);

	int isize, jsize, ksize;

private:
	Array3D<Material> grid;
};

inline bool getLineSegmentVoxelIntersection(Point3f p0, Point3f p1, double dx, CellMaterialGrid grid, Point3i* voxel) {
	double invdx = 1.0 / dx;
	p0 *= invdx;
	p1 *= invdx;

	int gx0idx = (int)std::floor(p0.x);
	int gy0idx = (int)std::floor(p0.y);
	int gz0idx = (int)std::floor(p0.z);

	int gx1idx = (int)std::floor(p1.x);
	int gy1idx = (int)std::floor(p1.y);
	int gz1idx = (int)std::floor(p1.z);

	int sx = gx1idx > gx0idx ? 1 : gx1idx < gx0idx ? -1 : 0;
	int sy = gy1idx > gy0idx ? 1 : gy1idx < gy0idx ? -1 : 0;
	int sz = gz1idx > gz0idx ? 1 : gz1idx < gz0idx ? -1 : 0;

	int gx = gx0idx;
	int gy = gy0idx;
	int gz = gz0idx;

	//Planes for each axis that we will next cross
	double gxp = gx0idx + (gx1idx > gx0idx ? 1 : 0);
	double gyp = gy0idx + (gy1idx > gy0idx ? 1 : 0);
	double gzp = gz0idx + (gz1idx > gz0idx ? 1 : 0);

	//Only used for multiplying up the error margins
	double vx = p1.x == p0.x ? 1 : p1.x - p0.x;
	double vy = p1.y == p0.y ? 1 : p1.y - p0.y;
	double vz = p1.z == p0.z ? 1 : p1.z - p0.z;

	//Error is normalized to vx * vy * vz so we only have to multiply up
	double vxvy = vx * vy;
	double vxvz = vx * vz;
	double vyvz = vy * vz;

	//Error from the next plane accumulators, scaled up by vx*vy*vz
	double errx = (gxp - p0.x) * vyvz;
	double erry = (gyp - p0.y) * vxvz;
	double errz = (gzp - p0.z) * vxvy;

	double derrx = sx * vyvz;
	double derry = sy * vxvz;
	double derrz = sz * vxvy;

	int gw = grid.isize;
	int gh = grid.jsize;
	int gd = grid.ksize;

	int maxiter = 1e6;
	int itercount = 0;
	for (;;) {
		if (isGridIndexInRange(gx, gy, gz, gw, gh, gd)) {
			if (grid.isCellSolid(gx, gy, gz)) {
				(*voxel).x = gx;
				(*voxel).y = gy;
				(*voxel).z = gz;
				return true;
			}
		}

		if (gx == gx1idx && gy == gy1idx && gz == gz1idx) {
			break;
		}

		//Which plane do we cross first?
		double xr = fabs(errx);
		double yr = fabs(erry);
		double zr = fabs(errz);

		if (sx != 0 && (sy == 0 || xr < yr) && (sz == 0 || xr < zr)) {
			gx += sx;
			errx += derrx;
		} else if (sy != 0 && (sz == 0 || yr < zr)) {
			gy += sy;
			erry += derry;
		} else if (sz != 0) {
			gz += sz;
			errz += derrz;
		}

		itercount++;
	}

	return false;
}

} // namespace foc

#endif // FOC_CELLMATERIALGRID_H
