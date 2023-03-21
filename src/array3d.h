#ifndef FOC_ARRAY3D_H
#define FOC_ARRAY3D_H

#include "foc.h"
#include "vecmath.h"

namespace foc {

template <typename T>
class Array3D {
public:
	Array3D() :
			width(0), height(0), depth(0), numElements(0) {
		initializeGrid();
	}

	Array3D(int i, int j, int k) :
			width(i), height(j), depth(k), numElements(i * j * k) {
		initializeGrid();
	}

	Array3D(int i, int j, int k, T fillval) :
			width(i), height(j), depth(k), numElements(i * j * k) {
		initializeGrid();
		fill(fillval);
	}

	Array3D operator=(const Array3D& arr) {
		delete[] grid;

		width = arr.width;
		height = arr.height;
		depth = arr.depth;
		numElements = arr.numElements;

		initializeGrid();

		T val;
		for (int k = 0; k < depth; k++) {
			for (int j = 0; j < height; j++) {
				for (int i = 0; i < width; i++) {
					val = arr.grid[getFlattenedIndex(i, j, k)];
					set(i, j, k, val);
				}
			}
		}

		return *this;
	}

	~Array3D() {
		delete[] grid;
	}

	T operator()(int i, int j, int k) {
		if (!isIndexInRange(i, j, k)) {
			std::printf("Error: index out of range.\n");
		}

		return grid[getFlattenedIndex(i, j, k)];
	}

	T operator()(int idx) {
		if (!(idx >= 0 && idx < numElements)) {
			std::printf("Error: index out of range.\n");
		}

		return grid[idx];
	}

	void fill(T value) {
		for (int idx = 0; idx < width * height * depth; idx++) {
			grid[idx] = value;
		}
	}

	void set(int i, int j, int k, T value) {
		if (!isIndexInRange(i, j, k)) {
			std::printf("Error: index out of range.\n");
		}

		grid[getFlattenedIndex(i, j, k)] = value;
	}

	void set(int idx, T value) {
		if (!(idx >= 0 && idx < numElements)) {
			std::printf("Error: index out of range.\n");
		}

		grid[idx] = value;
	}

	T* data() {
		return grid;
	}

	inline bool isIndexInRange(int i, int j, int k) {
		return i >= 0 && j >= 0 && k >= 0 && i < width && j < height && k < depth;
	}

	int width, height, depth;

private:
	void initializeGrid() {
		grid = new T[width * height * depth];
	}

	inline unsigned int getFlattenedIndex(int i, int j, int k) {
		return (unsigned int)i + (unsigned int)width * ((unsigned int)j + (unsigned int)height * (unsigned int)k);
	}

	T* grid;
	int numElements;
};

inline Point3f gridIndexToPosition(int i, int j, int k, double cellsize) {
	return Point3f(static_cast<double>(i) * cellsize, static_cast<double>(j) * cellsize, static_cast<double>(k) * cellsize);
}

inline Point3i positionToGridIndex(double x, double y, double z, double cellsize) {
	double invcs = 1.0 / cellsize;
	return Point3i(static_cast<int>(x * invcs), static_cast<int>(y * invcs), static_cast<int>(z * invcs));
}

inline Point3f gridIndexToCellCenter(int i, int j, int k, double cellsize) {
	double hcs = 0.5 * cellsize;
	return Point3f(static_cast<double>(i) * cellsize + hcs, static_cast<double>(j) * cellsize + hcs, static_cast<double>(k) * cellsize + hcs);
}

inline void getGridIndexBounds(Point3f pos, double r, double cellsize, int width, int height, int depth, Point3i* gmin, Point3i* gmax) {
	Point3i c = positionToGridIndex(pos.x, pos.y, pos.z, cellsize);
	Point3f cp = gridIndexToPosition(c.x, c.y, c.z, cellsize);
	Vector3f offset = pos - cp;
	double invcs = 1.0 / cellsize;

	int gimin = c.x - (int)fmax(0, std::ceil((r - offset.x) * invcs));
	int gjmin = c.y - (int)fmax(0, std::ceil((r - offset.y) * invcs));
	int gkmin = c.z - (int)fmax(0, std::ceil((r - offset.z) * invcs));
	int gimax = c.x + (int)fmax(0, std::ceil((r - cellsize + offset.x) * invcs));
	int gjmax = c.y + (int)fmax(0, std::ceil((r - cellsize + offset.y) * invcs));
	int gkmax = c.z + (int)fmax(0, std::ceil((r - cellsize + offset.z) * invcs));

	*gmin = Point3i((int)fmax(gimin, 0),
			(int)fmax(gjmin, 0),
			(int)fmax(gkmin, 0));
	*gmax = Point3i((int)fmin(gimax, width - 1),
			(int)fmin(gjmax, height - 1),
			(int)fmin(gkmax, depth - 1));
}

inline void getGridIndexBounds(Bounds3f bbox, double cellsize, int width, int height, int depth, Point3i* gmin, Point3i* gmax) {
	Vector3f offset = bbox.diagonal();
	*gmin = positionToGridIndex(bbox.p_min.x, bbox.p_min.y, bbox.p_min.z, cellsize);
	*gmax = positionToGridIndex(bbox.p_min.x + offset.x, bbox.p_min.y + offset.y, bbox.p_min.z + offset.z, cellsize);

	*gmin = Point3i((int)fmax((*gmin).x, 0),
			(int)fmax((*gmin).y, 0),
			(int)fmax((*gmin).z, 0));
	*gmax = Point3i((int)fmin((*gmax).x, width-1),
			(int)fmin((*gmax).y, height-1),
			(int)fmin((*gmax).z, depth-1));
}

inline void getGridIndexVertices(int i, int j, int k, Point3i v[8]) {
	v[0] = Point3i(i,     j,     k);
	v[1] = Point3i(i + 1, j,     k);
	v[2] = Point3i(i + 1, j,     k + 1);
	v[3] = Point3i(i,     j,     k + 1);
	v[4] = Point3i(i,     j + 1, k);
	v[5] = Point3i(i + 1, j + 1, k);
	v[6] = Point3i(i + 1, j + 1, k + 1);
	v[7] = Point3i(i,     j + 1, k + 1);
}

} // namespace foc

#endif // FOC_ARRAY3D_H
