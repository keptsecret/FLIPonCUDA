#ifndef FOC_TRIANGLEMESH_H
#define FOC_TRIANGLEMESH_H

#include "array3d.h"
#include "foc.h"

namespace foc {

class TriangleMesh {
public:
	TriangleMesh() {}

	int numVertices();
	int numFaces();
	void clear();

	void writeMeshToPLY(std::string filename);
	void writeMeshToBOBJ(std::string filename);

	void updateVertexTriangles();
	void updateVertexNormals();

	void smooth(double value, int iterations);
	void smooth(double value, int iterations, std::vector<int>& vertices);
	void translate(Vector3f offset);

	void removeMinimumTriangleCountPolyhedra(int count);

	std::vector<Point3f> vertices;
	std::vector<Point3f> vertexcolors; // r, g, b values in range [0.0, 1.0]
	std::vector<Vector3f> normals;
	std::vector<Point3i> indices;

private:
	void smoothTriangleMesh(double value, std::vector<bool> isVertexSmooth);

	int gridi = 0;
	int gridj = 0;
	int gridk = 0;
	double cellsize = 0;

	std::vector<std::vector<int>> vertexTriangles;
	std::vector<double> triangleAreas;

	Array3D<std::vector<int>> triGrid;
};

} // namespace foc

#endif // FOC_TRIANGLEMESH_H
