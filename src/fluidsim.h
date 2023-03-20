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

private:
	void initializeGrids(int i, int j, int k, int cellsize);

private:
	bool isInitialized = false;
	double gravity = -9.81;
	double density = 1000.0;
	double ratioFLIPPIC = 0.95;

	int isize, jsize, ksize;
	double dcell;

	std::vector<MarkerParticle> markerParticles;
	MACGrid macGrid;
	CellMaterialGrid materialGrid;
};

} // namespace foc

#endif // FOC_FLUIDSIM_H
