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

void TriangleMesh::writeMeshToPLY(std::string filename) {
	// Header format:
	/*
		ply
		format binary_little_endian 1.0
		element vertex FILL_IN_NUMBER_OF_VERTICES
		property float x
		property float y
		property float z
		element face FILL_IN_NUMBER_OF_FACES
		property list uchar int vertex_index
		end_header
	*/

	char header1[51] = { 'p', 'l', 'y', '\n',
		'f', 'o', 'r', 'm', 'a', 't', ' ', 'b', 'i', 'n', 'a', 'r', 'y', '_', 'l',
		'i', 't', 't', 'l', 'e', '_', 'e', 'n', 'd', 'i', 'a', 'n', ' ', '1', '.', '0', '\n',
		'e', 'l', 'e', 'm', 'e', 'n', 't', ' ', 'v', 'e', 'r', 't', 'e', 'x', ' ' };

	char header2[65] = { '\n', 'p', 'r', 'o', 'p', 'e', 'r', 't', 'y', ' ', 'f', 'l', 'o', 'a', 't', ' ', 'x', '\n',
		'p', 'r', 'o', 'p', 'e', 'r', 't', 'y', ' ', 'f', 'l', 'o', 'a', 't', ' ', 'y', '\n',
		'p', 'r', 'o', 'p', 'e', 'r', 't', 'y', ' ', 'f', 'l', 'o', 'a', 't', ' ', 'z', '\n',
		'e', 'l', 'e', 'm', 'e', 'n', 't', ' ', 'f', 'a', 'c', 'e', ' ' };

	char header2color[125] = { '\n', 'p', 'r', 'o', 'p', 'e', 'r', 't', 'y', ' ', 'f', 'l', 'o', 'a', 't', ' ', 'x', '\n',
		'p', 'r', 'o', 'p', 'e', 'r', 't', 'y', ' ', 'f', 'l', 'o', 'a', 't', ' ', 'y', '\n',
		'p', 'r', 'o', 'p', 'e', 'r', 't', 'y', ' ', 'f', 'l', 'o', 'a', 't', ' ', 'z', '\n',
		'p', 'r', 'o', 'p', 'e', 'r', 't', 'y', ' ', 'u', 'c', 'h', 'a', 'r', ' ', 'r', 'e', 'd', '\n',
		'p', 'r', 'o', 'p', 'e', 'r', 't', 'y', ' ', 'u', 'c', 'h', 'a', 'r', ' ', 'g', 'r', 'e', 'e', 'n', '\n',
		'p', 'r', 'o', 'p', 'e', 'r', 't', 'y', ' ', 'u', 'c', 'h', 'a', 'r', ' ', 'b', 'l', 'u', 'e', '\n',
		'e', 'l', 'e', 'm', 'e', 'n', 't', ' ', 'f', 'a', 'c', 'e', ' ' };

	char header3[49] = { '\n', 'p', 'r', 'o', 'p', 'e', 'r', 't', 'y', ' ', 'l', 'i', 's', 't', ' ',
		'u', 'c', 'h', 'a', 'r', ' ', 'i', 'n', 't', ' ',
		'v', 'e', 'r', 't', 'e', 'x', '_', 'i', 'n', 'd', 'e', 'x', '\n',
		'e', 'n', 'd', '_', 'h', 'e', 'a', 'd', 'e', 'r', '\n' };

	bool isColorEnabled = vertices.size() == vertexcolors.size();

	std::string vertstring = std::to_string(vertices.size());
	std::string facestring = std::to_string(indices.size());
	int vertdigits = (int)vertstring.length();
	int facedigits = (int)facestring.length();

	int offset = 0;
	int headersize;
	if (isColorEnabled) {
		headersize = 51 + vertdigits + 125 + facedigits + 49;
	} else {
		headersize = 51 + vertdigits + 65 + facedigits + 49;
	}

	int binsize;
	if (isColorEnabled) {
		binsize = headersize + 3 * (sizeof(float) * (int)vertices.size() + sizeof(unsigned char) * (int)vertices.size()) + (sizeof(unsigned char) + 3 * sizeof(int)) * (int)indices.size();
	} else {
		binsize = headersize + 3 * sizeof(float) * (int)vertices.size() + (sizeof(unsigned char) + 3 * sizeof(int)) * (int)indices.size();
	}
	char* bin = new char[binsize];

	memcpy(bin + offset, header1, 51);
	offset += 51;
	memcpy(bin + offset, vertstring.c_str(), vertdigits * sizeof(char));
	offset += vertdigits * sizeof(char);

	if (isColorEnabled) {
		memcpy(bin + offset, header2color, 125);
		offset += 125;
	} else {
		memcpy(bin + offset, header2, 65);
		offset += 65;
	}

	memcpy(bin + offset, facestring.c_str(), facedigits * sizeof(char));
	offset += facedigits * sizeof(char);
	memcpy(bin + offset, header3, 49);
	offset += 49;

	if (isColorEnabled) {
		float* vertdata = new float[3 * vertices.size()];
		Point3f v;
		for (unsigned int i = 0; i < vertices.size(); i++) {
			v = vertices[i];
			vertdata[3 * i] = v.x;
			vertdata[3 * i + 1] = v.y;
			vertdata[3 * i + 2] = v.z;
		}

		unsigned char* colordata = new unsigned char[3 * vertexcolors.size()];
		Point3f c;
		for (unsigned int i = 0; i < vertexcolors.size(); i++) {
			c = vertexcolors[i];
			colordata[3 * i] = (unsigned char)((c.x / 1.0) * 255.0);
			colordata[3 * i + 1] = (unsigned char)((c.y / 1.0) * 255.0);
			colordata[3 * i + 2] = (unsigned char)((c.z / 1.0) * 255.0);
		}

		int vertoffset = 0;
		int coloroffset = 0;
		int vertsize = 3 * sizeof(float);
		int colorsize = 3 * sizeof(unsigned char);
		for (unsigned int i = 0; i < vertices.size(); i++) {
			memcpy(bin + offset, vertdata + vertoffset, vertsize);
			offset += vertsize;
			vertoffset += 3;

			memcpy(bin + offset, colordata + coloroffset, colorsize);
			offset += colorsize;
			coloroffset += 3;
		}

		delete[] colordata;
		delete[] vertdata;
	} else {
		float* vertdata = new float[3 * vertices.size()];
		Point3f v;
		for (unsigned int i = 0; i < vertices.size(); i++) {
			v = vertices[i];
			vertdata[3 * i] = v.x;
			vertdata[3 * i + 1] = v.y;
			vertdata[3 * i + 2] = v.z;
		}
		memcpy(bin + offset, vertdata, 3 * sizeof(float) * vertices.size());
		offset += 3 * sizeof(float) * (int)vertices.size();
		delete[] vertdata;
	}

	int verts[3];
	for (unsigned int i = 0; i < indices.size(); i++) {
		Point3i t = indices[i];
		verts[0] = t[0];
		verts[1] = t[1];
		verts[2] = t[2];

		bin[offset] = 0x03;
		offset += sizeof(unsigned char);

		memcpy(bin + offset, verts, 3 * sizeof(int));
		offset += 3 * sizeof(int);
	}

	std::ofstream erasefile;
	erasefile.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc);
	erasefile.close();

	std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary);
	file.write(bin, binsize);
	file.close();

	delete[] bin;
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
		std::vector<int> triangles; // 14 is the maximum number of adjacent triangles to a vertex
		triangles.reserve(14);
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

	std::vector<Vector3f> faceNormals;
	faceNormals.reserve(indices.size());
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
	std::vector<int> verts;
	verts.reserve(vertices.size());
	for (int i = 0; i < vertices.size(); i++) {
		verts.push_back(i);
	}

	smooth(value, iterations, verts);
}

void TriangleMesh::smooth(double value, int iterations, std::vector<int>& verts) {
	std::vector<bool> isVertexSmooth(vertices.size(), false);
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
	std::vector<Point3f> newVertices;
	newVertices.reserve(vertices.size());

	for (int i = 0; i < vertices.size(); i++) {
		if (!isVertexSmooth[i]) {
			newVertices.push_back(vertices[i]);
			continue;
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
		Point3f nv = ov + (float)value * (avg - ov);
		newVertices.push_back(nv);
	}

	vertices = newVertices;
}

} // namespace foc