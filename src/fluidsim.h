#ifndef FOC_FLUIDSIM_H
#define FOC_FLUIDSIM_H

#include "foc.h"
#include "macgrid.h"
#include "cellmaterialgrid.h"
#include "markerparticle.h"

namespace foc {

class FluidSimulation {
public:
	FluidSimulation();

	FluidSimulation(int width, int height, int depth, int cellsize);

	~FluidSimulation() {}

	void initialize();
	void simulate(double fps, int numFrames);
	void update(double dt);

	void addFluidPoint(double x, double y, double z, double r);
	void addFluidPoint(Point3f p, double r);
	void addFluidCuboid(Point3f pmin, Point3f pmax);
	void addFluidCuboid(double x, double y, double z, double w, double h, double d);

private:
	void initializeGrids(int i, int j, int k, int cellsize);
	void initializeVectors(int i, int j, int k);

	void initializeSolidCells();
	void initializeFluid();
	void initializeMarkerParticles();

	void addMarkerParticlesToCell(Point3i idx);
	void addMarkerParticlesToCell(Point3i idx, Vector3f vel);

	void removeParticlesInSolidCells();
	void updateFluidCells();

	void computeVelocityField(Array3D<float>& field, Array3D<bool>& isValueSet, int dir);
	void advectVelocityFieldU();
	void advectVelocityFieldV();
	void advectVelocityFieldW();
	void advectVelocityField();

	double getNextTimeStep();
	double getMaxParticleSpeed();

	void stepSimulation(double dt);

private:
	struct FluidPoint {
		Point3f position;
		double radius = 0.0;

		FluidPoint() {}
		FluidPoint(Point3f pos, double r) : position(pos), radius(r) {}
	};

	struct FluidCuboid {
		Bounds3f bbox;

		FluidCuboid() {}
		FluidCuboid(Point3f pmin, Point3f pmax) : bbox(pmin, pmax) {}
	};

	bool isInitialized = false;
	double gravity = -9.81;
	double density = 1000.0;
	double ratioFLIPPIC = 0.95;
	FOC_CONST double CFLConditionNumber = 5.0;		// maximum number of cells a particle can move
	int randomSeed = 42;

	double markerParticleRadius = 0.0;
	int maxParticlesPerVelocityAdvection = 5e6;

	int isize, jsize, ksize;
	double dcell;

	bool isCurrentFrameFinished = false;
	int currentFrame = 0;
	double simulationTime = 0.0;

	std::vector<FluidPoint> fluidPoints;
	std::vector<FluidCuboid> fluidCuboids;

	std::vector<Point3i> fluidCellIndices;
	std::vector<MarkerParticle> markerParticles;
	MACGrid macGrid;
	CellMaterialGrid materialGrid;
};

} // namespace foc

#endif // FOC_FLUIDSIM_H
