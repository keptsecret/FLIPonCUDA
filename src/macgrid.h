#ifndef FOC_MACGRID_H
#define FOC_MACGRID_H

#include "foc.h"
#include "array3d.h"

namespace foc {

class MACGrid {
public:
	MACGrid();

	MACGrid(int width, int height, int depth, int cellsize);

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

private:
	void initializeVelocityGrids();

	int isize, jsize, ksize;
	double dcell;
	Array3D<float> u, v, w;
};

} // namespace foc

#endif // FOC_MACGRID_H
