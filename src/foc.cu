#include "foc.h"

#include "fluidsim.h"

int main() {
    foc::FluidSimulation fluidsim(64, 64, 64, 0.125);

	fluidsim.addFluidPoint(32, 32, 32, 5.0);
	fluidsim.addBodyForce(0, -9.81, 0);

	fluidsim.initialize();
	fluidsim.simulate(30.0, 30);

    return 0;
}
