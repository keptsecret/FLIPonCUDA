#include "particleadvector.h"

namespace foc {

ParticleAdvector::ParticleAdvector() {
}

void ParticleAdvector::tricubicInterpolate(std::vector<Vector3f>& particles, MACGrid* vfield, std::vector<Vector3f>& output) {
#ifndef FOC_BUILD_GPU
	output.resize(particles.size());

	for (int i = 0; i < particles.size(); i++) {
		output[i] = vfield->getVelocityAt(particles[i].x, particles[i].y, particles[i].z);
	}

	validateOutput(output);
#else
	// TODO: temporary code, change to implement gpu
	output.resize(particles.size());

	for (int i = 0; i < particles.size(); i++) {
		output[i] = vfield->getVelocityAt(particles[i].x, particles[i].y, particles[i].z);
	}

	validateOutput(output);
#endif
}

void ParticleAdvector::tricubicInterpolate(std::vector<Vector3f>& particles, MACGrid* vfield) {
	tricubicInterpolate(particles, vfield, particles);
}

void ParticleAdvector::advectParticles(std::vector<Point3f>& particles, MACGrid* vfield, double dt, std::vector<Point3f>& output) {
#ifndef FOC_BUILD_GPU
	output.clear();
	output.reserve(particles.size());
	for (int i = 0; i < particles.size(); i++) {
		output.push_back(RK4(particles[i], dt, vfield));
	}
#else
	output.clear();
	output.reserve(particles.size());

	std::vector<Vector3f> temppos;
	temppos.reserve(particles.size());
	for (int i = 0; i < particles.size(); i++) {
		output.push_back(particles[i]);
		temppos.push_back(Vector3f(particles[i]));
	}

	std::vector<Vector3f> tempdata;
	tempdata.reserve(particles.size());

	tricubicInterpolate(temppos, vfield, tempdata);

	float factor = dt / 6.0;
	for (int i = 0; i < tempdata.size(); i++) {
		output[i] += tempdata[i] * factor;
		tempdata[i] = temppos[i] + tempdata[i] * (0.5 * dt);
	}

	tricubicInterpolate(tempdata, vfield);

	for (int i = 0; i < tempdata.size(); i++) {
		output[i] += tempdata[i] * factor * 2.0;
		tempdata[i] = temppos[i] + tempdata[i] * (0.5 * dt);
	}

	tricubicInterpolate(tempdata, vfield);

	for (int i = 0; i < tempdata.size(); i++) {
		output[i] += tempdata[i] * factor * 2.0;
		tempdata[i] = temppos[i] + tempdata[i] * dt;
	}

	tricubicInterpolate(tempdata, vfield);

	for (int i = 0; i < tempdata.size(); i++) {
		output[i] += tempdata[i] * factor;
	}
#endif
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
	Vector3f k4 = vfield->getVelocityAt(p + k3 * dt);

	return p + (k1 + 2.0f * k2 + 2.0f * k3 + k4) * (dt / 6.0f);
}

} // namespace foc