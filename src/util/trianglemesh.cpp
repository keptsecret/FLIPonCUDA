#include "trianglemesh.h"
#include <fstream>

namespace foc {

int TriangleMesh::numVertices() {
	return vertices.size();
}

int TriangleMesh::numFaces() {
	return indices.size();
}

void TriangleMesh::clear() {
	vertices.clear();
	normals.clear();
	indices.clear();
	vertexcolors.clear();
	vertexTriangles.clear();
}

void TriangleMesh::writeMeshToBOBJ(std::string filename) {
	std::ofstream erasefile;
	erasefile.open(filename, std::ofstream::out | std::ofstream::trunc);
	erasefile.close();

	std::ofstream bobj(filename.c_str(), std::ios::out | std::ios::binary);

	int numVertices = (int)vertices.size();
	bobj.write((char*)&numVertices, sizeof(int));

	int binsize = 3 * numVertices * sizeof(float);
	bobj.write((char*)vertices.data(), binsize);

	int numTriangles = (int)indices.size();
	bobj.write((char*)&numTriangles, sizeof(int));

	binsize = 3 * numTriangles * sizeof(int);
	bobj.write((char*)indices.data(), binsize);

	bobj.close();
}

void TriangleMesh::updateVertexTriangles() {
	vertexTriangles.clear();
	vertexTriangles.reserve(vertices.size());

	for (int i = 0; i < vertices.size(); i++) {
		std::vector<int> triangles(14);	// 14 is the maximum number of adjacent triangles to a vertex
		vertexTriangles.push_back(triangles);
	}

	for (int i = 0; i < indices.size(); i++) {
		Point3i t = indices[i];
		vertexTriangles[t[0]].push_back(i);
		vertexTriangles[t[1]].push_back(i);
		vertexTriangles[t[2]].push_back(i);
	}
}

void TriangleMesh::updateVertexNormals() {
	normals.clear();
	updateVertexTriangles();

	std::vector<Vector3f> faceNormals(indices.size());
	for (int i = 0; i < indices.size(); i++) {
		Point3i t = indices[i];

		Vector3f v1 = vertices[t[1]] - vertices[t[0]];
		Vector3f v2 = vertices[t[2]] - vertices[t[0]];
		Vector3f normal = normalize(cross(v1, v2));

		faceNormals.push_back(normal);
	}

	for (int i = 0; i < vertexTriangles.size(); i++) {
		Vector3f n;
		for (int j = 0; j < vertexTriangles[i].size(); j++) {
			n += faceNormals[vertexTriangles[i][j]];
		}

		n = normalize(n / (float)vertexTriangles.size());
		normals.push_back(n);
	}
}

void TriangleMesh::smooth(double value, int iterations) {
	std::vector<int> verts(vertices.size());
	for (int i = 0; i < vertices.size(); i++) {
		verts.push_back(i);
	}

	smooth(value, iterations, verts);
}

void TriangleMesh::smooth(double value, int iterations, std::vector<int>& verts) {
	std::vector<bool> isVertexSmooth;
	isVertexSmooth.assign(vertices.size(), false);
	for (int i = 0; i < verts.size(); i++) {
		isVertexSmooth[verts[i]] = true;
	}

	vertexTriangles.clear();
	updateVertexTriangles();
	for (int i = 0; i < iterations; i++) {
		smoothTriangleMesh(value, isVertexSmooth);
	}
	vertexTriangles.clear();

	updateVertexNormals();
}

void TriangleMesh::translate(Vector3f offset) {
	for (int i = 0; i < vertices.size(); i++) {
		vertices[i] += offset;
	}
}

void TriangleMesh::removeMinimumTriangleCountPolyhedra(int count) {
	if (count <= 0) {
		return;
	}

	// TODO: might do something with this later
	//	std::vector<std::vector<int> > polyList;
	//	getPolyhedra(polyList);
	//
	//	std::vector<int> removalTriangles;
	//	for (unsigned int i = 0; i < polyList.size(); i++) {
	//		if ((int)polyList[i].size() <= count) {
	//			for (unsigned int j = 0; j < polyList[i].size(); j++) {
	//				removalTriangles.push_back(polyList[i][j]);
	//			}
	//		}
	//	}
	//
	//	if (removalTriangles.size() == 0) {
	//		return;
	//	}
	//
	//	removeTriangles(removalTriangles);
	//	removeExtraneousVertices();
}

void TriangleMesh::smoothTriangleMesh(double value, std::vector<bool> isVertexSmooth) {
	std::vector<Point3f> newVertices(vertices.size());

	for (int i = 0; i < vertices.size(); i++) {
		if (!isVertexSmooth[i]) {
			newVertices.push_back(vertices[i]);
		}

		int count = 0;
		Point3f avg;
		for (int j = 0; j < vertexTriangles[i].size(); j++) {
			Point3i t = indices[vertexTriangles[i][j]];
			if (t[0] != i) {
				avg += vertices[t[0]];
				count++;
			}
			if (t[1] != i) {
				avg += vertices[t[1]];
				count++;
			}
			if (t[2] != i) {
				avg += vertices[t[2]];
				count++;
			}
		}

		avg /= (float)count;
		Point3f ov = vertices[i];
		Point3f nv = ov + (float) value * (avg - ov);
		newVertices.push_back(nv);
	}

	vertices = newVertices;
}

} // namespace foc