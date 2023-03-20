#ifndef FOC_ARRAY3D_H
#define FOC_ARRAY3D_H

#include "foc.h"

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
		numElements = arr._numElements;

		initializeGrid();

		T val;
		for (int k = 0; k < depth; k++) {
			for (int j = 0; j < height; j++) {
				for (int i = 0; i < width; i++) {
					val = arr._grid[getFlattenedIndex(i, j, k)];
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

} // namespace foc

#endif // FOC_ARRAY3D_H
