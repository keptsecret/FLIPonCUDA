#include "fluidsim.h"

namespace foc {

FluidSimulation::FluidSimulation() {
}

FluidSimulation::FluidSimulation(int width, int height, int depth, int cellsize) :
		isize(width), jsize(height), ksize(depth), dcell(cellsize) {
	initializeGrids(width, height, depth, cellsize);
}

void FluidSimulation::initialize() {
	// TODO: initialize simulation here, e.g. mark fluid cells and create particles
}

void FluidSimulation::initializeGrids(int i, int j, int k, int cellsize) {
	macGrid = MACGrid(i, j, k, cellsize)
	materialGrid = CellMaterialGrid(i, j, k);
}

} // namespace foc