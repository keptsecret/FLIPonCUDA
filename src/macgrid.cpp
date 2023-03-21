#include "macgrid.h"

namespace foc {

MACGrid::MACGrid() :
		isize(10), jsize(10), ksize(10), dcell(0.1) {
	initializeVelocityGrids();
}

MACGrid::MACGrid(int width, int height, int depth, int cellsize) :
		isize(width), jsize(height), ksize(depth), dcell(cellsize) {
	initializeVelocityGrids();
}

void MACGrid::initializeVelocityGrids() {
	u = Array3D<float>(isize + 1, jsize, ksize, 0.0f);
	v = Array3D<float>(isize, jsize + 1, ksize, 0.0f);
	w = Array3D<float>(isize, jsize, ksize + 1, 0.0f);
}

void MACGrid::clearU() {
	u.fill(0.0f);
}

void MACGrid::clearV() {
	v.fill(0.0f);
}

void MACGrid::clearW() {
	w.fill(0.0f);
}

void MACGrid::clear() {
	clearU();
	clearV();
	clearW();
}

Array3D<float>* MACGrid::getArrayU() {
	return &u;
}

Array3D<float>* MACGrid::getArrayV() {
	return &v;
}

Array3D<float>* MACGrid::getArrayW() {
	return &w;
}

float* MACGrid::getRawArrayU() {
	return u.data();
}

float* MACGrid::getRawArrayV() {
	return v.data();
}

float* MACGrid::getRawArrayW() {
	return w.data();
}

float MACGrid::U(int i, int j, int k) {
	return u(i, j, k);
}

float MACGrid::V(int i, int j, int k) {
	return v(i, j, k);
}

float MACGrid::W(int i, int j, int k) {
	return w(i, j, k);
}

void MACGrid::setU(int i, int j, int k, double val) {
	u.set(i, j, k, val);
}

void MACGrid::setV(int i, int j, int k, double val) {
	v.set(i, j, k, val);
}

void MACGrid::setW(int i, int j, int k, double val) {
	w.set(i, j, k, val);
}

Vector3f MACGrid::getVelocityCell(int i, int j, int k) {
	double xa = 0.5 * (U(i + 1, j, k) + U(i, j, k));
	double ya = 0.5 * (V(i, j + 1, k) + V(i, j, k));
	double za = 0.5 * (W(i, j, k + 1) + W(i, j, k));

	return Vector3f(xa, ya, za);
}

float MACGrid::getVelocityMagCell(int i, int j, int k) {
	double magSq = getVelocityMagSqCell(i, j, k);
	if (magSq > 0.0) {
		return static_cast<float>(sqrt(magSq));
	}

	return 0.f;
}

float MACGrid::getVelocityMagSqCell(int i, int j, int k) {
	double xa = 0.5 * (U(i + 1, j, k) + U(i, j, k));
	double ya = 0.5 * (V(i, j + 1, k) + V(i, j, k));
	double za = 0.5 * (W(i, j, k + 1) + W(i, j, k));

	return static_cast<float>(xa * xa + ya * ya + za * za);
}

float MACGrid::getMaxVelocityMag() {
	double maxSq = 0.0;
	for (int k = 0; k < ksize; k++) {
		for (int j = 0; j < jsize; j++) {
			for (int i = 0; i < isize; i++) {
				double m = getVelocityMagSqCell(i, j, k);
				maxSq = fmax(maxSq, m);
			}
		}
	}

	double max = maxSq;
	if (maxSq > 0.0) {
		max = sqrt(maxSq);
	}

	return static_cast<float>(max);
}

Vector3f MACGrid::getVelocityFaceU(int i, int j, int k) {
	// Shift reference coordinate to the left. The formula used is for calculating
	// u(i+1/2, j, k). If we keep original (i,j,k) coordinate, then using the formula
	// would calculate u(i+3/2, j, k) instead. The same will be done for the V and W
	// faces, shifting back in the respective direction.

	i--;

	double vx = U(i + 1, j, k);
	double vy = 0.25 * (V(i, j, k) + V(i, j + 1, k) + V(i + 1, j, k) + V(i + 1, j + 1, k));
	double vz = 0.25 * (W(i, j, k) + W(i, j, k + 1) + W(i + 1, j, k) + W(i + 1, j, k + 1));

	return Vector3f(vx, vy, vz);
}

Vector3f MACGrid::getVelocityFaceV(int i, int j, int k) {
	j--;

	double vx = 0.25 * (U(i, j, k) + U(i + 1, j, k) + U(i, j + 1, k) + U(i + 1, j + 1, k));
	double vy = V(i, j + 1, k);
	double vz = 0.25 * (W(i, j, k) + W(i, j, k + 1) + W(i, j + 1, k) + W(i, j + 1, k + 1));

	return Vector3f(vx, vy, vz);
}

Vector3f MACGrid::getVelocityFaceW(int i, int j, int k) {
	k--;

	double vx = 0.25 * (U(i, j, k) + U(i + 1, j, k) + U(i, j, k + 1) + U(i + 1, j, k + 1));
	double vy = 0.25 * (V(i, j, k) + V(i, j + 1, k) + V(i, j, k + 1) + V(i, j + 1, k + 1));
	double vz = W(i, j, k + 1);

	return Vector3f(vx, vy, vz);
}

Vector3f MACGrid::getVelocityAt(double x, double y, double z) {
	// TODO: implement?
//	double xvel = interpolateU(x, y, z);
//	double yvel = interpolateV(x, y, z);
//	double zvel = interpolateW(x, y, z);

	return Vector3f();
}

Vector3f MACGrid::getVelocityAt(Point3f pos) {
	return getVelocityAt(pos.x, pos.y, pos.z);
}

} // namespace foc