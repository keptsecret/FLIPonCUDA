#include "fluidsim.h"

#include "scalarfield.h"
#include "pressuresolver.h"

namespace foc {

FluidSimulation::FluidSimulation() :
		isize(0), jsize(0), ksize(0), dcell(0) {
}

FluidSimulation::FluidSimulation(int width, int height, int depth, int cellsize) :
		isize(width), jsize(height), ksize(depth), dcell(cellsize) {
	initializeGrids(width, height, depth, cellsize);
	initializeVectors(width, height, depth);
}

void FluidSimulation::initialize() {
	// TODO: initialize simulation here, e.g. mark fluid cells and create particles
	if (!isInitialized) {
		srand(randomSeed);

		// initialize solid cells
		initializeSolidCells();

		// initialize fluid material cells
		initializeFluid();

		// initialize marker particles and radius
		initializeMarkerParticles();

		isInitialized = true;
	}
}

void FluidSimulation::simulate(double fps, int numFrames) {
	if (!isInitialized) {
		printf("Error: FluidSimluation has not been initialized. Aborting...");
		return;
	}

	double dt = 1.0 / fps;
	for (int i = 0; i < numFrames; i++) {
		update(dt);
	}
}

void FluidSimulation::update(double dt) {
	isCurrentFrameFinished = false;

	double timeleft = dt;
	while (timeleft > 0.0) {
		double timestep = getNextTimeStep();
		if (timeleft - timestep < 0.0) {
			timestep = timeleft;
		}
		timeleft -= timestep;

		stepSimulation(timestep);
	}

	currentFrame++;
	isCurrentFrameFinished = true;
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

void FluidSimulation::initializeVectors(int i, int j, int k) {
	fluidCellIndices = std::vector<Point3i>();
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
	std::vector<Point3i> fluidCells;
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

	// set fluid material grid and indices
	for (auto& mp : markerParticles) {
		Point3i idx = positionToGridIndex(mp.position.x, mp.position.y, mp.position.z, dcell);
		if (!materialGrid.isCellFluid(idx.x, idx.y, idx.z)) {
			materialGrid.set(idx.x, idx.y, idx.z, Material::fluid);
			fluidCellIndices.push_back(idx);
		}
	}
}

void FluidSimulation::initializeMarkerParticles() {
	double pvolume = dcell * dcell * dcell / 8.0;
	markerParticleRadius = pow((3 * pvolume) / (4 * foc::Pi), 1.0 / 3.0);
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

	double jitter = 0.1; // TODO: currently, constant jitter
	for (int i = 0; i < 8; i++) {
		Vector3f jit = Vector3f(randomDouble(-jitter, jitter),
				randomDouble(-jitter, jitter),
				randomDouble(-jitter, jitter));

		Point3f p = points[i] + jit;
		markerParticles.push_back(MarkerParticle(p, vel));
	}
}

void FluidSimulation::removeParticlesInSolidCells() {
	// removes marker particles in solid cells
	std::vector<bool> isRemoved(markerParticles.size());

	bool hasParticleInSolidCell = false;
	for (auto& mp : markerParticles) {
		Point3i idx = positionToGridIndex(mp.position.x, mp.position.y, mp.position.z, dcell);
		bool isInSolidCell = materialGrid.isCellSolid(idx.x, idx.y, idx.z);
		if (isInSolidCell) {
			hasParticleInSolidCell = true;
		}

		isRemoved.push_back(isInSolidCell);
	}

	// removes elements at positions true
	markerParticles.erase(std::remove_if(markerParticles.begin(), markerParticles.end(),
								  [&isRemoved, &markerParticles = markerParticles](auto const& i) {
									  return isRemoved.at(&i - markerParticles.data());
								  }),
			markerParticles.end());
	markerParticles.shrink_to_fit();
}

void FluidSimulation::updateFluidCells() {
	removeParticlesInSolidCells();

	for (int k = 1; k < materialGrid.ksize - 1; k++) {
		for (int j = 1; j < materialGrid.jsize - 1; j++) {
			for (int i = 1; i < materialGrid.isize - 1; i++) {
				if (materialGrid.isCellFluid(i, j, k)) {
					materialGrid.set(i, j, k, Material::air);
				}
			}
		}
	}

	fluidCellIndices.clear();

	// TODO: potential for parallelization
	for (const auto& mp : markerParticles) {
		Point3i idx = positionToGridIndex(mp.position.x, mp.position.y, mp.position.z, dcell);
		materialGrid.set(idx.x, idx.y, idx.z, Material::fluid);
	}

	int count = 0;
	for (int k = 0; k < materialGrid.ksize; k++) {
		for (int j = 0; j < materialGrid.jsize; j++) {
			for (int i = 0; i < materialGrid.isize; i++) {
				if (materialGrid.isCellFluid(i, j, k)) {
					count++;
				}
			}
		}
	}

	fluidCellIndices.reserve(count);
	for (int k = 0; k < materialGrid.ksize; k++) {
		for (int j = 0; j < materialGrid.jsize; j++) {
			for (int i = 0; i < materialGrid.isize; i++) {
				if (materialGrid.isCellFluid(i, j, k)) {
					fluidCellIndices.emplace_back(i, j, k);
				}
			}
		}
	}
}

void FluidSimulation::computeVelocityField(Array3D<float>& field, Array3D<bool>& isValueSet, int dir) {
	ScalarField velocityGrid(field.width, field.height, field.depth, dcell);
	velocityGrid.setPointRadius(dcell);
	velocityGrid.enableWeightField();

	Vector3f offset;
	switch (dir) {
		case 0: // U
			offset = Vector3f(0.0, 0.5 * dcell, 0.5 * dcell);
			break;
		case 1: // V
			offset = Vector3f(0.5 * dcell, 0.0, 0.5 * dcell);
			break;
		case 2: // U
			offset = Vector3f(0.5 * dcell, 0.5 * dcell, 0.0);
			break;
		default:
			return;
	}
	velocityGrid.setFieldOffset(offset);

	std::vector<Point3f> positions;
	std::vector<float> velocity; // velocity in direction
	positions.reserve(fmin(maxParticlesPerVelocityAdvection, markerParticles.size()));
	velocity.reserve(fmin(maxParticlesPerVelocityAdvection, markerParticles.size()));

	// TODO: potential for parallelization
#ifndef FOC_BUILD_GPU
	for (const auto& mp : markerParticles) {
		velocityGrid.addPointValue(mp.position, mp.velocity[dir]); // TODO: potentially check weights
	}
#else
	// for cuda acceleration here
#endif
	velocityGrid.applyWeightField();

	const auto scalarField = velocityGrid.getScalarFieldPointer();
	const auto weightField = velocityGrid.getWeightFieldPointer();
	for (int k = 0; k < field.depth; k++) {
		for (int j = 0; j < field.height; j++) {
			for (int i = 0; i < field.width; i++) {
				field.set(i, j, k, scalarField->get(i, j, k));

				if (weightField->get(i, j, k) > ShadowEpsilon) {
					isValueSet.set(i, j, k, true);
				}
			}
		}
	}
}

void FluidSimulation::advectVelocityFieldU() {
	macGrid.clearU();

	Array3D<float> uvel = Array3D<float>(isize + 1, jsize, ksize, 0.0f);
	Array3D<bool> isValueSet = Array3D<bool>(isize + 1, jsize, ksize, false);
	computeVelocityField(uvel, isValueSet, 0);

	std::vector<Point3i> extrapolationIndices((isize + 1) * jsize * ksize);
	for (int k = 0; k < uvel.depth; k++) {
		for (int j = 0; j < uvel.height; j++) {
			for (int i = 0; i < uvel.width; i++) {
				if (materialGrid.isFaceBorderingMaterialU(i, j, k, Material::fluid)) {
					if (!isValueSet(i, j, k)) {
						extrapolationIndices.push_back(Point3i(i, j, k));
					} else {
						macGrid.setU(i, j, k, uvel(i, j, k));
					}
				}
			}
		}
	}

	Point3i neighbors[26];
	for (const auto& g : extrapolationIndices) {
		getNeighborGridIndices26(g.x, g.y, g.z, neighbors);

		double avg = 0.0, weight = 0.0;
		for (int i = 0; i < 26; i++) {
			Point3i n = neighbors[i];
			if (uvel.isIndexInRange(n.x, n.y, n.z) && fabs(uvel(n.x, n.y, n.z)) > 0.0) {
				avg += uvel(n.x, n.y, n.z);
				weight += 1.0;
			}
		}

		if (weight > 0.0) {
			macGrid.setU(g.x, g.y, g.z, avg / weight);
		}
	}
}

void FluidSimulation::advectVelocityFieldV() {
	macGrid.clearV();

	Array3D<float> vvel = Array3D<float>(isize, jsize + 1, ksize, 0.0f);
	Array3D<bool> isValueSet = Array3D<bool>(isize, jsize + 1, ksize, false);
	computeVelocityField(vvel, isValueSet, 1);

	std::vector<Point3i> extrapolationIndices(isize * (jsize + 1) * ksize);
	for (int k = 0; k < vvel.depth; k++) {
		for (int j = 0; j < vvel.height; j++) {
			for (int i = 0; i < vvel.width; i++) {
				if (materialGrid.isFaceBorderingMaterialV(i, j, k, Material::fluid)) {
					if (!isValueSet(i, j, k)) {
						extrapolationIndices.push_back(Point3i(i, j, k));
					} else {
						macGrid.setV(i, j, k, vvel(i, j, k));
					}
				}
			}
		}
	}

	Point3i neighbors[26];
	for (const auto& g : extrapolationIndices) {
		getNeighborGridIndices26(g.x, g.y, g.z, neighbors);

		double avg = 0.0, weight = 0.0;
		for (int i = 0; i < 26; i++) {
			Point3i n = neighbors[i];
			if (vvel.isIndexInRange(n.x, n.y, n.z) && fabs(vvel(n.x, n.y, n.z)) > 0.0) {
				avg += vvel(n.x, n.y, n.z);
				weight += 1.0;
			}
		}

		if (weight > 0.0) {
			macGrid.setV(g.x, g.y, g.z, avg / weight);
		}
	}
}

void FluidSimulation::advectVelocityFieldW() {
	macGrid.clearW();

	Array3D<float> wvel = Array3D<float>(isize, jsize, ksize + 1, 0.0f);
	Array3D<bool> isValueSet = Array3D<bool>(isize, jsize, ksize + 1, false);
	computeVelocityField(wvel, isValueSet, 2);

	std::vector<Point3i> extrapolationIndices(isize * jsize * (ksize + 1));
	for (int k = 0; k < wvel.depth; k++) {
		for (int j = 0; j < wvel.height; j++) {
			for (int i = 0; i < wvel.width; i++) {
				if (materialGrid.isFaceBorderingMaterialW(i, j, k, Material::fluid)) {
					if (!isValueSet(i, j, k)) {
						extrapolationIndices.push_back(Point3i(i, j, k));
					} else {
						macGrid.setW(i, j, k, wvel(i, j, k));
					}
				}
			}
		}
	}

	Point3i neighbors[26];
	for (const auto& g : extrapolationIndices) {
		getNeighborGridIndices26(g.x, g.y, g.z, neighbors);

		double avg = 0.0, weight = 0.0;
		for (int i = 0; i < 26; i++) {
			Point3i n = neighbors[i];
			if (wvel.isIndexInRange(n.x, n.y, n.z) && fabs(wvel(n.x, n.y, n.z)) > 0.0) {
				avg += wvel(n.x, n.y, n.z);
				weight += 1.0;
			}
		}

		if (weight > 0.0) {
			macGrid.setW(g.x, g.y, g.z, avg / weight);
		}
	}
}

void FluidSimulation::updatePressureGrid(Array3D<float>& pressureGrid, double dt) {
	PressureSolverParameters params;
	params.cellsize = dcell;
	params.density = density;
	params.deltaTime = dt;
	params.fluidCells = &fluidCellIndices;
	params.materialGrid = &materialGrid;
	params.velocityField = &macGrid;

	VectorXd pressures(fluidCellIndices.size());
	PressureSolver solver;
	solver.solve(params, pressures);

	for (int i = 0; i < fluidCellIndices.size(); i++) {
		Point3i idx = fluidCellIndices[i];
		pressureGrid.set(idx.x, idx.y, idx.z, (float)pressures[i]);
	}
}

void FluidSimulation::applyPressureToVelocityField(Array3D<float>& pressureGrid, double dt) {
	// TODO: implement
}

void FluidSimulation::advectVelocityField() {
	advectVelocityFieldU();
	advectVelocityFieldV();
	advectVelocityFieldW();
}

double FluidSimulation::getNextTimeStep() {
	double maxv = getMaxParticleSpeed();
	double ts = CFLConditionNumber * dcell / maxv;
	return ts;
}

double FluidSimulation::getMaxParticleSpeed() {
	double maxvSq = 0.0;
	for (auto& mp : markerParticles) {
		double distSq = dot(mp.velocity, mp.velocity);
		if (distSq > maxvSq) {
			maxvSq = distSq;
		}
	}

	return sqrt(maxvSq);
}

void FluidSimulation::stepSimulation(double dt) {
	simulationTime += dt;

	// update grid with new marker particle positions and mark fluid cells
	updateFluidCells();

	// TODO: reconstruct level set (SDF) and meshes (usings OpenVDB)

	// advect velocity field with new marker particles
	advectVelocityField();

	// update and apply pressure
	Array3D<float> pressureGrid(isize, jsize, ksize, 0.0f);
	updatePressureGrid(pressureGrid, dt);
	applyPressureToVelocityField(pressureGrid, dt);

	// update (advect) marker particles
}

} // namespace foc