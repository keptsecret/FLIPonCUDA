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

	Array3D(const Array3D<T>& arr) {
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

		if (arr.bIsOutOfRangeValueSet) {
			OUT_OF_RANGE_VALUE = arr.OUT_OF_RANGE_VALUE;
			bIsOutOfRangeValueSet = true;
		}
	}

	Array3D<T>& operator=(const Array3D<T>& arr) {
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

		if (arr.bIsOutOfRangeValueSet) {
			OUT_OF_RANGE_VALUE = arr.OUT_OF_RANGE_VALUE;
			bIsOutOfRangeValueSet = true;
		}

		return *this;
	}

	~Array3D() {
		delete[] grid;
	}

	T operator()(int i, int j, int k) {
		if (!isIndexInRange(i, j, k)) {
			if (bIsOutOfRangeValueSet) {
				return OUT_OF_RANGE_VALUE;
			}
			std::printf("Error: index out of range.\n");
		}

		return grid[getFlattenedIndex(i, j, k)];
	}

	T operator()(int idx) {
		if (!(idx >= 0 && idx < numElements)) {
			if (bIsOutOfRangeValueSet) {
				return OUT_OF_RANGE_VALUE;
			}
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

	T get(int i, int j, int k) {
		if (!isIndexInRange(i, j, k)) {
			if (bIsOutOfRangeValueSet) {
				return OUT_OF_RANGE_VALUE;
			}
			std::printf("Error: index out of range.\n");
		}

		return grid[getFlattenedIndex(i, j, k)];
	}

	T get(int idx) {
		if (!(idx >= 0 && idx < numElements)) {
			if (bIsOutOfRangeValueSet) {
				return OUT_OF_RANGE_VALUE;
			}
			std::printf("Error: index out of range.\n");
		}

		return grid[idx];
	}

	T* getPointer(int i, int j, int k) {
		if (!isIndexInRange(i, j, k)) {
			if (bIsOutOfRangeValueSet) {
				return &OUT_OF_RANGE_VALUE;
			}
			std::printf("Error: index out of range.\n");
		}

		return &grid[getFlattenedIndex(i, j, k)];
	}

	T* getPointer(int idx) {
		if (!(idx >= 0 && idx < numElements)) {
			if (bIsOutOfRangeValueSet) {
				return &OUT_OF_RANGE_VALUE;
			}
			std::printf("Error: index out of range.\n");
		}

		return &grid[idx];
	}

	T* data() {
		return grid;
	}

	void setOutOfRangeValue() {
		bIsOutOfRangeValueSet = false;
	}
	void setOutOfRangeValue(T val) {
		OUT_OF_RANGE_VALUE = val;
		bIsOutOfRangeValueSet = true;
	}

	bool isOutOfRangeValueSet() {
		return bIsOutOfRangeValueSet;
	}

	T getOutOfRangeValue() {
		return OUT_OF_RANGE_VALUE;
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

	bool bIsOutOfRangeValueSet = false;
	T OUT_OF_RANGE_VALUE;

	T* grid;
	int numElements;
};

inline bool isGridIndexInRange(Point3i g, Point3i gmax) {
	return g.x >= 0 && g.y >= 0 && g.z >= 0 && g.x <= gmax.x && g.y <= gmax.y && g.z <= gmax.z;
}

inline bool isGridIndexInRange(Point3i g, int imax, int jmax, int kmax) {
	return g.x >= 0 && g.y >= 0 && g.z >= 0 && g.x <= imax && g.y <= jmax && g.z <= kmax;
}

inline bool isGridIndexInRange(int i, int j, int k, int imax, int jmax, int kmax) {
	return i >= 0 && j >= 0 && k >= 0 && i <= imax && j <= jmax && k <= kmax;
}

inline bool isPositionInGrid(double x, double y, double z, double cellsize, int i, int j, int k) {
	return x >= 0 && y >= 0 && z >= 0 && x < cellsize * i && y < cellsize * j && z < cellsize * k;
}

inline Point3f gridIndexToPosition(int i, int j, int k, double cellsize) {
	return Point3f(static_cast<double>(i) * cellsize, static_cast<double>(j) * cellsize, static_cast<double>(k) * cellsize);
}

inline Point3i positionToGridIndex(double x, double y, double z, double cellsize) {
	double invcs = 1.0 / cellsize;
	return Point3i(static_cast<int>(std::floor(x * invcs)), static_cast<int>(std::floor(y * invcs)), static_cast<int>(std::floor(z * invcs)));
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
	*gmax = Point3i((int)fmin((*gmax).x, width - 1),
			(int)fmin((*gmax).y, height - 1),
			(int)fmin((*gmax).z, depth - 1));
}

inline void getGridIndexVertices(int i, int j, int k, Point3i v[8]) {
	v[0] = Point3i(i, j, k);
	v[1] = Point3i(i + 1, j, k);
	v[2] = Point3i(i + 1, j, k + 1);
	v[3] = Point3i(i, j, k + 1);
	v[4] = Point3i(i, j + 1, k);
	v[5] = Point3i(i + 1, j + 1, k);
	v[6] = Point3i(i + 1, j + 1, k + 1);
	v[7] = Point3i(i, j + 1, k + 1);
}

inline void getNeighborGridIndices6(int i, int j, int k, Point3i nb[6]) {
	nb[0] = Point3i(i - 1, j, k);
	nb[1] = Point3i(i + 1, j, k);
	nb[2] = Point3i(i, j - 1, k);
	nb[3] = Point3i(i, j + 1, k);
	nb[4] = Point3i(i, j, k - 1);
	nb[5] = Point3i(i, j, k + 1);
}

inline void getNeighborGridIndices26(int i, int j, int k, Point3i nb[26]) {
	int idx = 0;
	for (int nk = k - 1; nk <= k + 1; nk++) {
		for (int nj = j - 1; nj <= j + 1; nj++) {
			for (int ni = i - 1; ni <= i + 1; ni++) {
				if (!(ni == i && nj == j && nk == k)) {
					nb[idx] = Point3i(ni, nj, nk);
					idx++;
				}
			}
		}
	}
}

} // namespace foc

#endif // FOC_ARRAY3D_H
