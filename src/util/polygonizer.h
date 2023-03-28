#ifndef FOC_POLYGONIZER_H
#define FOC_POLYGONIZER_H

#include "array3d.h"
#include "foc.h"
#include "scalarfield.h"
#include "trianglemesh.h"

namespace foc {

class Polygonizer {
public:
	Polygonizer() {}

	Polygonizer(ScalarField* field);

	void setSurfaceCellMask(Array3D<bool>* mask);
	TriangleMesh polygonizeSurface();

private:
	struct EdgeGrid {
		Array3D<int> U; // store index to vertex
		Array3D<int> V;
		Array3D<int> W;

		EdgeGrid() :
				U(Array3D<int>(0, 0, 0)),
				V(Array3D<int>(0, 0, 0)),
				W(Array3D<int>(0, 0, 0)) {}

		EdgeGrid(int i, int j, int k) :
				U(Array3D<int>(i, j + 1, k + 1, -1)),
				V(Array3D<int>(i + 1, j, k + 1, -1)),
				W(Array3D<int>(i + 1, j + 1, k, -1)) {}
	};

	bool isCellOnSurface(Point3i idx);
	void findSurfaceCells(std::vector<Point3i>& surfaceCells);

	// Marching cubes to create mesh
	int calculateCubeIndex(Point3i cell);
	Point3f vertexInterpolate(Point3f p1, Point3f p2, double v1, double v2);
	void calculateVertexList(Point3i cell, int cubeIndex, EdgeGrid& edges, std::vector<Point3f>& meshVertices, int vertexList[12]);
	void polygonizeCell(Point3i cell, EdgeGrid& edges, TriangleMesh& mesh);
	void calculateSurfaceTriangles(std::vector<Point3i>& surfaceCells, TriangleMesh& mesh);

private:
	static const int edgeTable[256];
	static const int triTable[256][16];

	int isize = 0;
	int jsize = 0;
	int ksize = 0;
	double cellsize = 0.0;

	double surfaceThreshold = 0.5;

	ScalarField* scalarField;
	bool isScalarFieldSet = false;

	Array3D<bool>* surfaceCellMask;
	bool isSurfaceCellMaskSet = false;
};

} // namespace foc

#endif // FOC_POLYGONIZER_H
