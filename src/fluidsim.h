#ifndef FOC_FLUIDSIM_H
#define FOC_FLUIDSIM_H

#include "cellmaterialgrid.h"
#include "foc.h"
#include "macgrid.h"
#include "markerparticle.h"

namespace foc {

class FluidSimulation {
public:
	FluidSimulation();

	FluidSimulation(int width, int height, int depth, double cellsize);

	~FluidSimulation() {}

	void initialize();
	void simulate(double fps, int numFrames);
	void update(double dt);

	void addFluidPoint(double x, double y, double z, double r);
	void addFluidPoint(Point3f p, double r);
	void addFluidCuboid(Point3f pmin, Point3f pmax);
	void addFluidCuboid(double x, double y, double z, double w, double h, double d);

	// constant body force
	void addBodyForce(double fx, double fy, double fz);
	void addBodyForce(Vector3f f);
	void resetBodyForce();

	Vector3f getConstantBodyForce();

private:
	void initializeGrids(int i, int j, int k, int cellsize);
	void initializeVectors(int i, int j, int k);

	void initializeSolidCells();
	void initializeFluid();
	void initializeMarkerParticles();

	void addMarkerParticlesToCell(Point3i idx);
	void addMarkerParticlesToCell(Point3i idx, Vector3f vel);

	// update grid with new marker particle positions and mark fluid cells
	void removeParticlesInSolidCells();
	void updateFluidCells();

	// advect velocity field with new marker particles
	void computeVelocityField(Array3D<float>& field, Array3D<bool>& isValueSet, int dir);
	void advectVelocityFieldU();
	void advectVelocityFieldV();
	void advectVelocityFieldW();
	void advectVelocityField();

	// apply body forces, e.g. gravity
	void applyConstantBodyForces(double dt);
	void applyBodyForcesToVelocityField(double dt);

	// update and apply pressure
	void updatePressureGrid(Array3D<float>& pressureGrid, double dt);
	void applyPressureToFaceU(int i, int j, int k, Array3D<float>& pressureGrid, MACGrid& updatedMACGrid, double deltaTime);
	void applyPressureToFaceV(int i, int j, int k, Array3D<float>& pressureGrid, MACGrid& updatedMACGrid, double deltaTime);
	void applyPressureToFaceW(int i, int j, int k, Array3D<float>& pressureGrid, MACGrid& updatedMACGrid, double deltaTime);
	void applyPressureToVelocityField(Array3D<float>& pressureGrid, double dt);

	void extrapolateFluidVelocities(MACGrid& velocityGrid);

	// update (advect) marker particles
	void updateMarkerParticleVelocities();


	double getNextTimeStep();
	double getMaxParticleSpeed();

	void stepSimulation(double dt);

private:
	struct FluidPoint {
		Point3f position;
		double radius = 0.0;

		FluidPoint() {}
		FluidPoint(Point3f pos, double r) :
				position(pos), radius(r) {}
	};

	struct FluidCuboid {
		Bounds3f bbox;

		FluidCuboid() {}
		FluidCuboid(Point3f pmin, Point3f pmax) :
				bbox(pmin, pmax) {}
	};

	bool isInitialized = false;
	double gravity = -9.81;
	double density = 1000.0;
	FOC_CONST double CFLConditionNumber = 5.0; // maximum number of cells a particle can move
	int randomSeed = 42;

	double markerParticleRadius = 0.0;
	int maxParticlesPerVelocityAdvection = 5e6;
	double ratioFLIPPIC = 0.95;
	int maxParticlesPerFLIPPICUpdate = 10e6;

	int isize, jsize, ksize;
	double dcell;

	bool isCurrentFrameFinished = false;
	int currentFrame = 0;
	double simulationTime = 0.0;

	std::vector<FluidPoint> fluidPoints;
	std::vector<FluidCuboid> fluidCuboids;

	std::vector<Vector3f> constantBodyForces;

	std::vector<Point3i> fluidCellIndices;
	std::vector<MarkerParticle> markerParticles;
	MACGrid macGrid;
	CellMaterialGrid materialGrid;
};

} // namespace foc

#endif // FOC_FLUIDSIM_H
