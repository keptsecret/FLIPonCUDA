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

	void addFluidPoint(double x, double y, double z, double r);
	void addFluidPoint(Point3f p, double r);
	void addFluidCuboid(Point3f pmin, Point3f pmax);
	void addFluidCuboid(double x, double y, double z, double w, double h, double d);

private:
	void initializeGrids(int i, int j, int k, int cellsize);

	void initializeSolidCells();
	void initializeFluid();

	void addMarkerParticlesToCell(Point3i idx);
	void addMarkerParticlesToCell(Point3i idx, Vector3f vel);

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

	double markerParticleRadius = 0.0;

	int isize, jsize, ksize;
	double dcell;

	int currentFrame;

	std::vector<FluidPoint> fluidPoints;
	std::vector<FluidCuboid> fluidCuboids;
	std::vector<Point3i> fluidCells;

	std::vector<MarkerParticle> markerParticles;
	MACGrid macGrid;
	CellMaterialGrid materialGrid;
};

} // namespace foc

#endif // FOC_FLUIDSIM_H
