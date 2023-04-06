#ifndef FOC_ARRAY3DVIEW_H
#define FOC_ARRAY3DVIEW_H

#include "array3d.h"

namespace foc {

template <typename T>
class Array3DView {
public:
	Array3DView() {
		setDimensions(0, 0, 0);
		setOffset(0, 0, 0);
		setArray3D(&dummyGrid);
	}

	Array3DView(Array3D<T>* grid) {
		setDimensions(0, 0, 0);
		setOffset(0, 0, 0);
		setArray3D(grid);
	}

	Array3DView(int isize, int jsize, int ksize, Array3D<T>* grid) {
		setDimensions(isize, jsize, ksize);
		setOffset(0, 0, 0);
		setArray3D(grid);
	}

	Array3DView(int isize, int jsize, int ksize, int offi, int offj, int offk, Array3D<T>* grid) {
		setDimensions(isize, jsize, ksize);
		setOffset(offi, offj, offk);
		setArray3D(grid);
	}

	Array3DView(const Array3DView& arr) {
		width = arr.width;
		height = arr.height;
		depth = arr.depth;

		ioffset = arr.ioffset;
		joffset = arr.joffset;
		koffset = arr.koffset;

		dummyGrid = Array3D<T>();

		if (arr.parent == &arr.dummyGrid) {
			parent = &dummyGrid;
		} else {
			parent = arr.parent;
		}
	}

	Array3DView& operator=(const Array3DView& arr) {
		width = arr.width;
		height = arr.height;
		depth = arr.depth;

		ioffset = arr.ioffset;
		joffset = arr.joffset;
		koffset = arr.koffset;

		dummyGrid = Array3D<T>();

		if (arr.parent == &arr.dummyGrid) {
			parent = &dummyGrid;
		} else {
			parent = arr.parent;
		}

		return *this;
	}

	void setDimensions(int isize, int jsize, int ksize) {
		width = isize;
		height = jsize;
		depth = ksize;
	}

	void getDimensions(int* isize, int* jsize, int* ksize) {
		*isize = width;
		*ksize = height;
		*jsize = depth;
	}

	void setOffset(int offi, int offj, int offk) {
		ioffset = offi;
		joffset = offj;
		koffset = offk;
	}

	void setArray3D(Array3D<T>* grid) {
		parent = grid;
	}

	Array3D<T>* getArray3D() {
		return parent;
	}

	T get(int i, int j, int k) {
		if (!isIndexInView(i, j, k)) {
			std::printf("Error: index out of view range.\n");
		}

		Point3i parentIdx = viewToParentIndex(i, j, k);
		if (!parent->isIndexInRange(parentIdx.x, parentIdx.y, parentIdx.z)) {
			if (parent->isOutOfRangeValueSet()) {
				return parent->getOutOfRangeValue();
			}
			std::printf("Error: index out of range.\n");
		}

		return parent->get(parentIdx.x, parentIdx.y, parentIdx.z);
	}

	T operator()(int i, int j, int k) {
		return get(i, j, k);
	}

	void set(int i, int j, int k, T value) {
		if (!isIndexInView(i, j, k)) {
			std::printf("Error: index out of view range.\n");
		}

		Point3i parentIdx = viewToParentIndex(i, j, k);
		if (parent->isIndexInRange(parentIdx.x, parentIdx.y, parentIdx.z)) {
			parent->set(parentIdx.x, parentIdx.y, parentIdx.z, value);
		}
	}

	inline bool isIndexInView(int i, int j, int k) {
		return i >= 0 && j >= 0 && k >= 0 && i < width && j < height && k < depth;
	}

	inline bool isIndexInParent(int i, int j, int k) {
		Point3i pidx = viewToParentIndex(i, j, k);
		return parent->isIndexInRange(pidx.x, pidx.y, pidx.z);
	}

	int width, height, depth;

private:
	inline Point3i viewToParentIndex(int i, int j, int k) {
		return Point3i (i + ioffset, j + joffset, k + koffset);
	}

	int ioffset = 0;
	int joffset = 0;
	int koffset = 0;

	Array3D<T>* parent;
	Array3D<T> dummyGrid;
};

} // namespace foc

#endif // FOC_ARRAY3DVIEW_H
