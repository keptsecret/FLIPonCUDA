#ifndef FOC_SCALARFIELD_H
#define FOC_SCALARFIELD_H

#include "array3d.h"
#include "cellmaterialgrid.h"
#include "foc.h"

namespace foc {

class ScalarField {
public:
	ScalarField();

	ScalarField(int i, int j, int k, int cellsize);

	void setPointRadius(double r);
	void setFieldOffset(Vector3f offset);

	void enableWeightField();
	void applyWeightField();

	void addPoint(Point3f p, double r);
	void addPoint(Point3f p);
	void addCuboid(Point3f pmin, Point3f pmax);

	void addPointValue(Point3f p, double val);

	void setScalarFieldValue(int i, int j, int k, double value);
	void addScalarFieldValue(int i, int j, int k, double value);
	double getScalarFieldValue(int i, int j, int k);

	void getGridDimensions(int* i, int* j, int* k) {
		*i = isize;
		*j = jsize;
		*k = ksize;
	}

	double getCellSize() { return cellsize; }

	Array3D<float>* getScalarFieldPointer();
	Array3D<float>* getWeightFieldPointer();

	void setMaterialGrid(CellMaterialGrid& materialGrid);

	double getSurfaceThreshold() { return surfaceThreshold; }

	double getScalarFieldValueAtCellCenter(int i, int j, int k);

private:
	double evaluateTricubicFieldFunctionForRadiusSquared(double rsq);

	int isize = 0;
	int jsize = 0;
	int ksize = 0;
	double cellsize = 0.0;

	double radius = 0.0;
	double invRadius = 1.0;
	double coef1 = 0.0;
	double coef2 = 0.0;
	double coef3 = 0.0;

	double surfaceThreshold = 0.5;
	double maxScalarFieldThreshold = 0.0;
	bool isMaxScalarFieldThresholdSet = false;

	Array3D<float> field;
	Array3D<bool> isVertexSolid;
	Array3D<float> weightField;
	Array3D<bool> isVertexSet;

	bool isWeightFieldEnabled = false;

	Vector3f gridOffset;
};

} // namespace foc

#endif // FOC_SCALARFIELD_H
