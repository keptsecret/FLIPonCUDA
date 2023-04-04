#ifndef FOC_PARTICLEADVECTOR_H
#define FOC_PARTICLEADVECTOR_H

#include "foc.h"
#include "macgrid.h"
#include "util/array3dview.h"

namespace foc {

class ParticleAdvector {
public:
	ParticleAdvector();

	void tricubicInterpolate(std::vector<Vector3f>& particles, MACGrid* vfield, std::vector<Vector3f>& output);
	void tricubicInterpolate(std::vector<Vector3f>& particles, MACGrid* vfield);

	void advectParticles(std::vector<Point3f>& particles, MACGrid* vfield, double dt, std::vector<Point3f>& output);

private:
	void validateOutput(std::vector<Vector3f>& output);

	Point3f RK4(Point3f p, double dt, MACGrid* vfield);

	struct ParticleBatch {
		std::vector<Vector3f> particles;	// particle positions
		std::vector<int> refs;
	};

	struct BatchData {
		std::vector<Vector3f>::iterator particlesBegin;
		std::vector<Vector3f>::iterator particlesEnd;
		std::vector<int>::iterator refsBegin;
		std::vector<int>::iterator refsEnd;

		Array3DView<float> ufieldView;
		Array3DView<float> vfieldView;
		Array3DView<float> wfieldView;

		Point3i batchOffset;
		Point3i indexOffset;
		Vector3f positionOffset;
	};

	int getMaxBatchesPerComputation();
	void getParticleBatchGrid(double bw, double bh, double bd, std::vector<Vector3f>& particles, Array3D<ParticleBatch>& grid);
	void getBatchData(MACGrid* vfield, Array3D<ParticleBatch>& particleGrid, std::vector<BatchData>& batchData);

	void tricubicInterpolateBatch(std::vector<BatchData>& batches, std::vector<Vector3f>& output);


private:
	int isize = 0;
	int jsize = 0;
	int ksize = 0;
	double cellsize = 0.0;

	int batchDim = 5;
	int maxThreadsPerBlock = 512;	// block size, CUDA's limit is 1024
	int maxBatchesPerComputation = 15000; // grid size
};

} // namespace foc

#endif // FOC_PARTICLEADVECTOR_H
