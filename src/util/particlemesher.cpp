#include "particlemesher.h"
#include "polygonizer.h"

namespace foc {

ParticleMesher::ParticleMesher(int i, int j, int k, double dx) :
		isize(i), jsize(j), ksize(k), cellsize(dx) {
}

void ParticleMesher::setSubdivisionLevel(int n) {
	subdivisionLevel = n;
}

void ParticleMesher::setNumPolygonizationSlices(int n) {
	numPolygonizationSlices = n > isize ? isize : n;
}

TriangleMesh ParticleMesher::particleToMesh(std::vector<MarkerParticle>& particles, CellMaterialGrid& materialGrid, double radius) {
	particleRadius = radius;
	if (numPolygonizationSlices == 1) {
		return polygonizeAll(particles, materialGrid);
	}

	return polygonizeSlices(particles, materialGrid);
}

TriangleMesh ParticleMesher::polygonizeAll(std::vector<MarkerParticle>& particles, CellMaterialGrid& materialGrid) {
	int subd = subdivisionLevel;
	int width = isize * subd;
	int height = jsize * subd;
	int depth = ksize * subd;
	double dx = cellsize / (double)subd;

	ScalarField field(width + 1, height + 1, depth + 1, dx);

	field.setMaterialGrid(materialGrid);
	field.setPointRadius(particleRadius);
	addPointsToScalarField(particles, field);

	Polygonizer polygonizer(&field);

	return polygonizer.polygonizeSurface();
}

TriangleMesh ParticleMesher::polygonizeSlices(std::vector<MarkerParticle>& particles, CellMaterialGrid& materialGrid) {
	// TODO: unimplemented
	return TriangleMesh();
}

void ParticleMesher::addPointsToScalarField(std::vector<MarkerParticle>& particles, ScalarField& field) {
	for (int i = 0; i < particles.size(); i++) {
		field.addPoint(particles[i].position);
	}
}

} // namespace foc