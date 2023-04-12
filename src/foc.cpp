#include "foc.h"

#include "fluidsim.h"

int main() {

	// Fluid sphere
	foc::FluidSimulation fluidsim(64, 64, 64, 0.125);

	double width, height, depth;
	fluidsim.getSimulationDimensions(&width, &height, &depth);
	fluidsim.addFluidPoint(width / 2, height / 2, depth / 2, 6.0);
	fluidsim.addBodyForce(0, -25, 0);

	fluidsim.initialize();
	fluidsim.simulate(30.0, 60);

	// Fluid wave with obstacle
	/*
	foc::FluidSimulation fluidsim(128, 64, 64, 0.125);

	double width, height, depth;
	fluidsim.getSimulationDimensions(&width, &height, &depth);
	fluidsim.addFluidCuboid(0, 0, 0, depth * 0.75, height * 0.5, depth);
	fluidsim.addBodyForce(0, -25, 0);

	std::vector<foc::Point3i> solidBlock;
	for (int k = 16; k < 48; k++) {
		for (int j = 0; j < 16; j++) {
			for (int i = 84; i < 100; i++) {
				solidBlock.push_back(foc::Point3i(i, j, k));
			}
		}
	}
	fluidsim.addSolidCells(solidBlock);

	fluidsim.initialize();
	fluidsim.simulate(30.0, 120);
	 */

    return 0;
}
