#include "foc.h"

#include "fluidsim.h"

int main() {
    foc::FluidSimulation fluidsim(32, 32, 32, 0.25);

	double width, height, depth;
	fluidsim.getSimulationDimensions(&width, &height, &depth);
	fluidsim.addFluidPoint(width / 2, height / 2, depth / 2, 6.0);
	fluidsim.addBodyForce(0, -25, 0);

	fluidsim.initialize();
	fluidsim.simulate(30.0, 60);

    return 0;
}
