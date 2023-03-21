#include "fluidsim.h"

#include "scalarfield.h"

namespace foc {

FluidSimulation::FluidSimulation() {
}

FluidSimulation::FluidSimulation(int width, int height, int depth, int cellsize) :
		isize(width), jsize(height), ksize(depth), dcell(cellsize) {
	initializeGrids(width, height, depth, cellsize);
}

void FluidSimulation::initialize() {
	// TODO: initialize simulation here, e.g. mark fluid cells and create particles
	if (!isInitialized) {
		// initialize solid cells
		initializeSolidCells();

		// initialize fluid material cells
		initializeFluid();

		// initialize marker particles and radius

		isInitialized = true;
	}
}

void FluidSimulation::addFluidPoint(double x, double y, double z, double r) {
	addFluidPoint(Point3f(x, y, z), r);
}

void FluidSimulation::addFluidPoint(Point3f p, double r) {
	fluidPoints.emplace_back(p, r);
}

void FluidSimulation::addFluidCuboid(Point3f pmin, Point3f pmax) {
	fluidCuboids.emplace_back(pmin, pmax);
}

void FluidSimulation::addFluidCuboid(double x, double y, double z, double w, double h, double d) {
	addFluidCuboid(Point3f(x, y, z), Point3f(x + w, y + h, z + d));
}

void FluidSimulation::initializeGrids(int i, int j, int k, int cellsize) {
	macGrid = MACGrid(i, j, k, cellsize);
	materialGrid = CellMaterialGrid(i, j, k);
}

void FluidSimulation::initializeSolidCells() {
	// sets border cells to be solid
	for (int j = 0; j < jsize; j++) {
		for (int i = 0; i < isize; i++) {
			materialGrid.set(i, j, 0, Material::solid);
			materialGrid.set(i, j, ksize - 1, Material::solid);
		}
	}

	for (int k = 0; k < ksize; k++) {
		for (int i = 0; i < isize; i++) {
			materialGrid.set(i, 0, k, Material::solid);
			materialGrid.set(i, jsize - 1, k, Material::solid);
		}
	}

	for (int k = 0; k < ksize; k++) {
		for (int j = 0; j < jsize; j++) {
			materialGrid.set(0, j, k, Material::solid);
			materialGrid.set(isize - 1, j, k, Material::solid);
		}
	}
}

void FluidSimulation::initializeFluid() {
	// turn fluid points and cuboids into a scalar field
	ScalarField field = ScalarField(isize + 1, jsize + 1, ksize + 1, dcell);
	for (auto& fp : fluidPoints) {
		field.addPoint(fp.position, fp.radius);
	}
	for (auto& fc : fluidCuboids) {
		field.addCuboid(fc.bbox.p_min, fc.bbox.p_max);
	}

	// get fluid cells from scalar field
	field.setMaterialGrid(materialGrid);
	double threshold = field.getSurfaceThreshold();

	Vector3f c;
	for (int k = 0; k < ksize; k++) {
		for (int j = 0; j < jsize; j++) {
			for (int i = 0; i < isize; i++) {
				if (!materialGrid.isCellSolid(i, j, k) && field.getScalarFieldValueAtCellCenter(i, j, k) > threshold) {
					fluidCells.push_back(Point3i(i, j, k));
				}
			}
		}
	}

	// add marker particles to cell (TODO: currently disregards any partial cells)
	markerParticles.reserve(8 * fluidCells.size());
	for (int i = 0; i < fluidCells.size(); i++) {
		addMarkerParticlesToCell(fluidCells[i]);
	}
}

void FluidSimulation::addMarkerParticlesToCell(Point3i idx) {
	addMarkerParticlesToCell(idx, Vector3f());
}

void FluidSimulation::addMarkerParticlesToCell(Point3i idx, Vector3f vel) {
	double q = 0.25 * dcell;
	Point3f c = gridIndexToCellCenter(idx.x, idx.y, idx.z, dcell);

	Point3f points[] = {
		Point3f(c.x - q, c.y - q, c.z - q),
		Point3f(c.x + q, c.y - q, c.z - q),
		Point3f(c.x + q, c.y - q, c.z + q),
		Point3f(c.x - q, c.y - q, c.z + q),
		Point3f(c.x - q, c.y + q, c.z - q),
		Point3f(c.x + q, c.y + q, c.z - q),
		Point3f(c.x + q, c.y + q, c.z + q),
		Point3f(c.x - q, c.y + q, c.z + q)
	};

	double jitter = 0.1; 	// TODO: currently, constant jitter
	for (int idx = 0; idx < 8; idx++) {
		Vector3f jit = Vector3f(randomDouble(-jitter, jitter),
				randomDouble(-jitter, jitter),
				randomDouble(-jitter, jitter));

		Point3f p = points[idx] + jit;
		markerParticles.push_back(MarkerParticle(p, vel));
	}
}

} // namespace foc