#include "particleadvector.h"

namespace foc {

ParticleAdvector::ParticleAdvector() {
}

void ParticleAdvector::tricubicInterpolate(std::vector<Point3f>& particles, MACGrid* vfield, std::vector<Vector3f>& output) {
#ifndef FOC_BUILD_GPU
	output.resize(particles.size());

	for (int i = 0; i < particles.size(); i++) {
		output[i] = vfield->getVelocityAt(particles[i]);
	}

	validateOutput(output);
#endif
	return;
}

void ParticleAdvector::advectParticles(std::vector<Point3f>& particles, MACGrid* vfield, double dt, std::vector<Point3f>& output) {
	output.clear();
	output.reserve(particles.size());
	for (int i = 0; i < particles.size(); i++) {
		output.push_back(RK4(particles[i], dt, vfield));
	}
}

void ParticleAdvector::validateOutput(std::vector<Vector3f>& output) {
	for (auto& v : output) {
		if (std::isinf(v.x) || std::isnan(v.x) ||
				std::isinf(v.y) || std::isnan(v.y) ||
				std::isinf(v.z) || std::isnan(v.z)) {
			v = Vector3f();
		}
	}
}

Point3f ParticleAdvector::RK4(Point3f p, double dt, MACGrid* vfield) {
	Vector3f k1 = vfield->getVelocityAt(p);
	Vector3f k2 = vfield->getVelocityAt(p + k1 * 0.5 * dt);
	Vector3f k3 = vfield->getVelocityAt(p + k2 * 0.5 * dt);
	Vector3f k4 = vfield->getVelocityAt(p + k3 * 0.5 * dt);

	return p + (k1 + 2.0f * k2 + 2.0f * k3 + k4) * (dt / 6.0f);
}

} // namespace foc