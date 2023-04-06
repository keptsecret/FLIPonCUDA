#ifndef FOC_SCALARFIELD_H
#define FOC_SCALARFIELD_H

#include "array3d.h"
#include "array3dview.h"
#include "cellmaterialgrid.h"
#include "foc.h"

namespace foc {

class ScalarField {
public:
	ScalarField();

	ScalarField(int i, int j, int k, double cellsize);

	void setPointRadius(double r);
	void setFieldOffset(Vector3f offset);

	void enableWeightField();
	void applyWeightField();

	void addPoint(Point3f p, double r);
	void addPoint(Point3f p);
	void addCuboid(Point3f pmin, Point3f pmax);

	void addPointValue(Point3f p, double val);

	// calls GPU functions
	void addPointValues(std::vector<Vector3f>& points, std::vector<float>& values);
	void addPointValues(std::vector<Vector3f>& points, std::vector<float>& values, double radius, Vector3f offset, double dx);

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

	struct PointData {
		PointData() {}
		PointData(Vector3f p, float v) :
				position(p), value(v) {}

		Vector3f position;
		float value;
	};

	struct BatchData {
		std::vector<PointData> particles;

		Array3DView<float> fieldView;
		Array3DView<float> weightfieldView;

		Point3i batchOffset;
		Point3i indexOffset;
		Vector3f positionOffset;

		float minFieldValue = 0.0;
	};

	struct BatchDataRef {
		Point3i batchDataIndex;

		std::vector<PointData>::iterator particlesBegin;
		std::vector<PointData>::iterator particlesEnd;
	};

	void initializePointData(std::vector<Vector3f>& points, std::vector<float>& values, std::vector<PointData>& pd);
	void initializeBatchDataGrid(std::vector<PointData>& points, Array3D<BatchData>& grid);
	void initializeBatchRefs(Array3D<BatchData>& grid, std::vector<BatchDataRef>& queue);

	int getMaxBatchesPerPointDataComputation();
	void updateBatchMinimumValues(Array3D<BatchData>& grid);
	void getNextBatches(std::vector<BatchDataRef>& queue, Array3D<BatchData>& grid, std::vector<BatchDataRef>& batches, int n);

	void computePointValueScalarField(std::vector<BatchDataRef>& batches, Array3D<BatchData>& batchDataGrid);
	void computePointValueScalarWeightField(std::vector<BatchDataRef>& batches, Array3D<BatchData>& batchDataGrid);

	int getMaxNumParticlesInBatch(std::vector<BatchDataRef>& batches);
	void fillPointValueDataBuffer(std::vector<BatchDataRef>& batches, Array3D<BatchData>& grid, int numParticles, std::vector<float>& buffer);
	void fillScalarFieldDataBuffer(std::vector<BatchDataRef>& batches, std::vector<float>& buffer);
	void fillScalarWeightFieldDataBuffer(std::vector<BatchDataRef>& batches, std::vector<float>& buffer);
	void fillBatchOffsetDataBuffer(std::vector<BatchDataRef>& batches, std::vector<Point3i>& buffer);

private:
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

	int batchDim = 8; // cbrt(maxThreadsPerBlock)
	int maxThreadsPerBlock = 512; // block size, CUDA's limit is 1024
	int maxBatchesPerComputation = 15000; // grid size
	int maxParticlesPerBatch = 1000;

	Vector3f gridOffset;
};

} // namespace foc

#endif // FOC_SCALARFIELD_H
