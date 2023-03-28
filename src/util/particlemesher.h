#ifndef FOC_PARTICLEMESHER_H
#define FOC_PARTICLEMESHER_H

#include "foc.h"
#include "markerparticle.h"
#include "scalarfield.h"
#include "trianglemesh.h"

namespace foc {

class ParticleMesher {
public:
	ParticleMesher() {}

	ParticleMesher(int i, int j, int k, double dx);

	void setSubdivisionLevel(int n);
	void setNumPolygonizationSlices(int n);

	TriangleMesh particleToMesh(std::vector<MarkerParticle>& particles, CellMaterialGrid& materialGrid, double radius);

private:
	TriangleMesh polygonizeAll(std::vector<MarkerParticle>& particles, CellMaterialGrid& materialGrid);
	TriangleMesh polygonizeSlices(std::vector<MarkerParticle>& particles, CellMaterialGrid& materialGrid);

	void addPointsToScalarField(std::vector<MarkerParticle>& particles, ScalarField& field);

private:
	int isize = 0;
	int jsize = 0;
	int ksize = 0;
	double cellsize = 0.0;

	int subdivisionLevel = 1;
	int numPolygonizationSlices = 1;

	double particleRadius = 0.0;
	double maxScalarFieldValueThreshold = 1.0;

	Array3D<float> scalarFieldSeamData;

	int maxParticlesPerScalarFieldAddition = 5e6;

	bool isPreviewMesherEnabled = false;
	int pisize = 0;
	int pjsize = 0;
	int pksize = 0;
	double pcellsize = 0.0;
	ScalarField pfield;
};

} // namespace foc

#endif // FOC_PARTICLEMESHER_H
