#ifndef FOC_PARTICLEADVECTOR_H
#define FOC_PARTICLEADVECTOR_H

#include "foc.h"
#include "macgrid.h"

namespace foc {

class ParticleAdvector {
public:
	ParticleAdvector();

	//bool initialize();

	void tricubicInterpolate(std::vector<Vector3f>& particles, MACGrid* vfield, std::vector<Vector3f>& output);
	void tricubicInterpolate(std::vector<Vector3f>& particles, MACGrid* vfield);

	void advectParticles(std::vector<Point3f>& particles, MACGrid* vfield, double dt, std::vector<Point3f>& output);

private:
	void validateOutput(std::vector<Vector3f>& output);

	Point3f RK4(Point3f p, double dt, MACGrid* vfield);
};

} // namespace foc

#endif // FOC_PARTICLEADVECTOR_H
