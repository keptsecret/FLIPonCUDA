#include "scalarfield.h"

namespace foc {

ScalarField::ScalarField() {
}

ScalarField::ScalarField(int i, int j, int k, int cellsize) :
		isize(i), jsize(j), ksize(k), cellsize(cellsize), field(i, j, k, false), isVertexSolid(i, j, k), isVertexSet(i, j, k, false) {
}

void ScalarField::setPointRadius(double r) {
	radius = r;
	invRadius = 1 / r;
	coef1 = (4.0 / 9.0) * (1.0 / (r * r * r * r * r * r));
	coef2 = (17.0 / 9.0) * (1.0 / (r * r * r * r));
	coef3 = (22.0 / 9.0) * (1.0 / (r * r));
}

void ScalarField::setFieldOffset(Vector3f offset) {
	gridOffset = offset;
}

void ScalarField::enableWeightField() {
	if (isWeightFieldEnabled) {
		return;
	}

	weightField = Array3D<float>(isize, jsize, ksize, 0.0f);
	isWeightFieldEnabled = true;
}

void ScalarField::applyWeightField() {
	if (!isWeightFieldEnabled) {
		return;
	}

	for (int k = 0; k < ksize; k++) {
		for (int j = 0; j < jsize; j++) {
			for (int i = 0; i < isize; i++) {
				float weight = weightField(i, j, k);
				if (weight > 0.0) {
					float v = field(i, j, k) / weight;
					setScalarFieldValue(i, j, k, v);
				}
			}
		}
	}
}

void ScalarField::addPoint(Point3f p, double r) {
	setPointRadius(r);
	addPoint(p);
}

void ScalarField::addPoint(Point3f p) {
	p -= gridOffset;

	Point3i gmin, gmax;
	getGridIndexBounds(p, radius, cellsize, isize, jsize, ksize, &gmin, &gmax);

	Point3f gpos;
	Vector3f v;
	double rsq = radius * radius;
	double distsq;
	double weight;
	for (int k = gmin.z; k <= gmax.z; k++) {
		for (int j = gmin.y; j <= gmax.y; j++) {
			for (int i = gmin.x; i <= gmax.x; i++) {
				if (isMaxScalarFieldThresholdSet && field(i, j, k) > maxScalarFieldThreshold) {
					continue;
				}

				gpos = gridIndexToPosition(i, j, k, cellsize);
				v = gpos - p;
				distsq = dot(v, v);
				if (distsq < rsq) {
					weight = evaluateTricubicFieldFunctionForRadiusSquared(distsq);
					addScalarFieldValue(i, j, k, weight);

					if (isWeightFieldEnabled) {
						weightField.set(i, j, k, weightField(i, j, k) + (float)weight);
					}
				}
			}
		}
	}
}

void ScalarField::addCuboid(Point3f pmin, Point3f pmax) {
	pmin -= gridOffset;

	Point3i gmin, gmax;
	Bounds3f bbox = Bounds3f(pmin, pmax);
	getGridIndexBounds(bbox, cellsize, isize, jsize, ksize, &gmin, &gmax);


	Point3f gpos;
	for (int k = gmin.z; k <= gmax.z; k++) {
		for (int j = gmin.y; j <= gmax.y; j++) {
			for (int i = gmin.x; i <= gmax.x; i++) {
				if (isMaxScalarFieldThresholdSet && field(i, j, k) > maxScalarFieldThreshold) {
					continue;
				}

				gpos = gridIndexToPosition(i, j, k, cellsize);
				if (bInside(gpos, bbox)) {
					addScalarFieldValue(i, j, k, surfaceThreshold + foc::ShadowEpsilon);

					if (isWeightFieldEnabled) {
						weightField.set(i, j, k, (float)(weightField(i, j, k) + surfaceThreshold + foc::ShadowEpsilon));
					}
				}
			}
		}
	}
}

void ScalarField::addPointValue(Point3f p, double val) {
	p -= gridOffset;

	Point3i gmin, gmax;
	getGridIndexBounds(p, radius, cellsize, isize, jsize, ksize, &gmin, &gmax);

	Point3f gpos;
	Vector3f v;
	double rsq = radius * radius;
	double distsq;
	double weight;
	for (int k = gmin.z; k <= gmax.z; k++) {
		for (int j = gmin.y; j <= gmax.y; j++) {
			for (int i = gmin.x; i <= gmax.x; i++) {
				if (isMaxScalarFieldThresholdSet && field(i, j, k) > maxScalarFieldThreshold) {
					continue;
				}

				gpos = gridIndexToPosition(i, j, k, cellsize);
				v = gpos - p;
				distsq = dot(v, v);
				if (distsq < rsq) {
					weight = evaluateTricubicFieldFunctionForRadiusSquared(distsq);
					addScalarFieldValue(i, j, k, weight * val);

					if (isWeightFieldEnabled) {
						weightField.set(i, j, k, weightField(i, j, k) + (float)weight);
					}
				}
			}
		}
	}
}

void ScalarField::setScalarFieldValue(int i, int j, int k, double value) {
	field.set(i, j, k, value);
	isVertexSet.set(i, j, k, true);
}

void ScalarField::addScalarFieldValue(int i, int j, int k, double value) {
	field.set(i, j, k, field(i, j, k) + value);
	isVertexSet.set(i, j, k, true);
}

double ScalarField::getScalarFieldValue(int i, int j, int k) {
	double val = field(i, j, k);
	if (isVertexSolid(i, j, k) && val > surfaceThreshold) {
		val = surfaceThreshold;
	}

	return val;
}

Array3D<float>* ScalarField::getScalarFieldPointer() {
	return &field;
}

Array3D<float>* ScalarField::getWeightFieldPointer() {
	return &weightField;
}

void ScalarField::setMaterialGrid(CellMaterialGrid& materialGrid) {
	Point3i vertices[8];
	for (int k = 0; k < ksize - 1; k++) {
		for (int j = 0; j < jsize - 1; j++) {
			for (int i = 0; i < isize - 1; i++) {
				if (materialGrid.isCellSolid(i, j, k)) {
					getGridIndexVertices(i, j, k, vertices);
					for (int idx = 0; idx < 8; idx++) {
						isVertexSolid.set(vertices[idx].x, vertices[idx].y, vertices[idx].z, true);
					}
				}
			}
		}
	}
}

double ScalarField::getScalarFieldValueAtCellCenter(int i, int j, int k) {
	double sum = 0.0;
	sum += getScalarFieldValue(i, j, k);
	sum += getScalarFieldValue(i + 1, j, k);
	sum += getScalarFieldValue(i, j + 1, k);
	sum += getScalarFieldValue(i + 1, j + 1, k);
	sum += getScalarFieldValue(i, j, k + 1);
	sum += getScalarFieldValue(i + 1, j, k + 1);
	sum += getScalarFieldValue(i, j + 1, k + 1);
	sum += getScalarFieldValue(i + 1, j + 1, k + 1);

	return 0.125 * sum;
}

double ScalarField::evaluateTricubicFieldFunctionForRadiusSquared(double rsq) {
	return 1.0 - coef1 * rsq * rsq * rsq + coef2 * rsq * rsq - coef3 * rsq;
}

} // namespace foc