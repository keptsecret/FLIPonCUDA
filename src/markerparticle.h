#ifndef FOC_MARKERPARTICLE_H
#define FOC_MARKERPARTICLE_H

#include "foc.h"
#include "vecmath.h"

namespace foc {

class MarkerParticle {
public:
	MarkerParticle() {}

	MarkerParticle(Point3f pos) :
			position(pos) {}

	MarkerParticle(Point3f pos, Vector3f vel) :
			position(pos), velocity(vel) {}

	MarkerParticle(double x, double y, double z) :
			position(0.0, 0.0, 0.0) {}

	Point3f position;
	Vector3f velocity;
};

} // namespace foc

#endif // FOC_MARKERPARTICLE_H
