#ifndef FOC_FLUIDSIM_H
#define FOC_FLUIDSIM_H

#include "cellmaterialgrid.h"
#include "foc.h"
#include "macgrid.h"
#include "markerparticle.h"
#include "particleadvector.h"
#include "util/trianglemesh.h"

namespace foc {

class FluidSimulation {
public:
	FluidSimulation();

	FluidSimulation(int width, int height, int depth, double cellsize);

	~FluidSimulation() {}

	void getSimulationDimensions(double* width, double* height, double* depth);

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
	void initializeGrids(int i, int j, int k, double cellsize);
	void initializeVectors(int i, int j, int k);

	void initializeSolidCells();
	void initializeFluid();
	void initializeMarkerParticles();

	void addMarkerParticlesToCell(Point3i idx);
	void addMarkerParticlesToCell(Point3i idx, Vector3f vel);

	// update grid with new marker particle positions and mark fluid cells
	void removeParticlesInSolidCells();
	void updateFluidCells();

	// reconstruct level set (SDF) and meshes
	void polygonizeOutputSurface(TriangleMesh& mesh);
	bool isVertexNearSolid(Point3f p, double eps);
	void smoothSurfaceMesh(TriangleMesh& mesh);
	void reconstructFluidSurfaceMesh();
	void exportMeshToFile(TriangleMesh& mesh, std::string filename);

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
	void updateMarkerParticleVelocitiesSubset(int start, int end);
	void updateMarkerParticleVelocities();
	void advanceMarkerParticlesSubset(int start, int end, double dt);
	void removeMarkerParticles();
	Point3f resolveParticleSolidCollision(Point3f oldpos, Point3f newpos);
	void advanceMarkerParticles(double dt);

	double getNextTimeStep();
	double getMaxParticleSpeed();

	void stepSimulation(double dt);

private:
	struct FluidPoint {
		Point3f position;
		double radius = 0.0;

		FluidPoint() {}
		FluidPoint(Point3f pos, double r) :
				position(pos), radius(r) {
		}
	};

	struct FluidCuboid {
		Bounds3f bbox;

		FluidCuboid() {}
		FluidCuboid(Point3f pmin, Point3f pmax) :
				bbox(pmin, pmax) {}
	};

	bool isInitialized = false;
	double density = 20.0;
	double CFLConditionNumber = 5.0; // maximum number of cells a particle can move
	int randomSeed = 42;

	double markerParticleRadius = 0.0;
	double markerParticleScale = 3.0;
	int maxParticlesPerVelocityAdvection = 5e6;
	double ratioFLIPPIC = 0.95;
	int maxParticlesPerFLIPPICUpdate = 10e6;
	int maxMarkerParticlesPerCell = 100;

	int minimumSurfacePolyhedronTriangleCount = 0;
	double surfaceSmoothingValue = 0.5;
	int surfaceSmoothingIterations = 2;
	Vector3f domainOffset;

	int isize, jsize, ksize;
	double cellsize;

	bool isCurrentFrameFinished = false;
	int currentFrame = 0;
	double simulationTime = 0.0;

	std::vector<FluidPoint> fluidPoints;
	std::vector<FluidCuboid> fluidCuboids;

	std::vector<Vector3f> constantBodyForces;

	ParticleAdvector particleAdvector;

	std::vector<Point3i> fluidCellIndices;
	std::vector<MarkerParticle> markerParticles;
	MACGrid macGrid;
	MACGrid prevMACGrid;
	CellMaterialGrid materialGrid;

	// for logging times
	std::vector<std::chrono::duration<double>> times;
};

} // namespace foc

#endif // FOC_FLUIDSIM_H
