#include "fluidsim.h"

#include "pressuresolver.h"
#include "util/particlemesher.h"
#include "util/scalarfield.h"
#include <chrono>

namespace foc {

FluidSimulation::FluidSimulation() :
		isize(0), jsize(0), ksize(0), dcell(0) {
}

FluidSimulation::FluidSimulation(int width, int height, int depth, double cellsize) :
		isize(width), jsize(height), ksize(depth), dcell(cellsize) {
	initializeGrids(isize, jsize, ksize, dcell);
	initializeVectors(isize, jsize, ksize);
}

void FluidSimulation::getSimulationDimensions(double* width, double* height, double* depth) {
	*width = (double)isize * dcell;
	*height = (double)jsize * dcell;
	*depth = (double)ksize * dcell;
}

void FluidSimulation::initialize() {
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
	auto start = std::chrono::steady_clock::now();

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

	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "Frame " << currentFrame << ":: elapsed time: " << elapsed_seconds.count() << "s\n";
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

void FluidSimulation::addBodyForce(double fx, double fy, double fz) {
	addBodyForce(Vector3f(fx, fy, fz));
}

void FluidSimulation::addBodyForce(Vector3f f) {
	constantBodyForces.push_back(f);
}

void FluidSimulation::resetBodyForce() {
	constantBodyForces.clear();
}

Vector3f FluidSimulation::getConstantBodyForce() {
	Vector3f total;
	for (const auto& bf : constantBodyForces) {
		total += bf;
	}
	return total;
}

void FluidSimulation::initializeGrids(int i, int j, int k, double cellsize) {
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
	ScalarField field(isize + 1, jsize + 1, ksize + 1, dcell);
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
	std::vector<bool> isRemoved;
	isRemoved.reserve(markerParticles.size());

	bool hasParticleInSolidCell = false;
	for (auto& mp : markerParticles) {
		Point3i idx = positionToGridIndex(mp.position.x, mp.position.y, mp.position.z, dcell);
		bool isInSolidCell = materialGrid.isCellSolid(idx.x, idx.y, idx.z);
		if (isInSolidCell) {
			hasParticleInSolidCell = true;
		}

		isRemoved.push_back(isInSolidCell);
	}

	if (hasParticleInSolidCell) {
		// removes elements at positions true
		markerParticles.erase(std::remove_if(markerParticles.begin(), markerParticles.end(),
									  [&isRemoved, &markerParticles = markerParticles](auto const& i) {
										  return isRemoved.at(&i - markerParticles.data());
									  }),
				markerParticles.end());
		markerParticles.shrink_to_fit();
	}
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

void FluidSimulation::polygonizeOutputSurface(TriangleMesh& mesh) {
	double radius = markerParticleRadius * markerParticleScale;

	ParticleMesher mesher(isize, jsize, ksize, dcell);
	mesh = mesher.particleToMesh(markerParticles, materialGrid, radius);
}

bool FluidSimulation::isVertexNearSolid(Point3f p, double eps) {
	Point3i gridIndex = positionToGridIndex(p.x, p.y, p.z, dcell);
	if (materialGrid.isCellSolid(gridIndex.x, gridIndex.y, gridIndex.z)) {
		return true;
	}

	// is p near the solid boundary?
	Point3f e(eps, eps, eps);
	if (gridIndex.x == 1 || gridIndex.x == isize - 2 ||
			gridIndex.y == 1 || gridIndex.y == jsize - 2 ||
			gridIndex.z == 1 || gridIndex.z == ksize - 2) {
		Point3f min = gridIndexToPosition(1, 1, 1, dcell);
		Point3f max = gridIndexToPosition(isize - 2, jsize - 2, ksize - 2, dcell);
		Bounds3f bbox(min + e, max + Vector3f(dcell, dcell, dcell) - Vector3f(e));
		if (!bInside(p, bbox)) {
			return true;
		}
	}

	// is p near a solid cell?
	Point3f gp = gridIndexToPosition(gridIndex.x, gridIndex.y, gridIndex.z, dcell);
	Bounds3f bbox(gp + e, gp + Vector3f(dcell, dcell, dcell) - Vector3f(e));
	if (bInside(p, bbox)) {
		return false;
	}

	Vector3f ex(eps, 0.0, 0.0);
	Vector3f ey(0.0, eps, 0.0);
	Vector3f ez(0.0, 0.0, eps);

	Point3f points[26]{
		p - ex, p + ex, p - ey, p + ey, p - ez, p + ez, p - ex - ey, p - ex + ey,
		p + ex - ey, p + ex + ey, p - ex - ez, p - ex + ez, p + ex - ez, p + ex + ez,
		p - ey - ez, p - ey + ez, p + ey - ez, p + ey + ez, p - ex - ey - ez,
		p - ex - ey + ez, p - ex + ey - ez, p - ex + ey + ez, p + ex - ey - ez,
		p + ex - ey + ez, p + ex + ey - ez, p + ex + ey + ez
	};

	for (int idx = 0; idx < 26; idx++) {
		Point3i gidx = positionToGridIndex(points[idx].x, points[idx].y, points[idx].z, dcell);
		if (materialGrid.isCellSolid(gidx.x, gidx.y, gidx.z)) {
			return true;
		}
	}

	return false;
}

void FluidSimulation::smoothSurfaceMesh(TriangleMesh& mesh) {
	std::vector<int> smoothVertices;
	double eps = 0.02 * dcell;
	for (int i = 0; i < mesh.vertices.size(); i++) {
		Point3f v = mesh.vertices[i];
		if (!isVertexNearSolid(v, eps)) {
			smoothVertices.push_back(i);
		}
	}

	mesh.smooth(surfaceSmoothingValue, surfaceSmoothingIterations, smoothVertices);
}

void FluidSimulation::exportMeshToFile(TriangleMesh& mesh, std::string filename) {
	mesh.writeMeshToPLY(filename);
}

void FluidSimulation::reconstructFluidSurfaceMesh() {
	TriangleMesh isoMesh;
	polygonizeOutputSurface(isoMesh);
	isoMesh.removeMinimumTriangleCountPolyhedra(minimumSurfacePolyhedronTriangleCount);

	smoothSurfaceMesh(isoMesh);
	isoMesh.translate(domainOffset);

	std::string framestr = std::to_string(currentFrame);
	framestr.insert(framestr.begin(), 6 - framestr.size(), '0');
	std::string bakedir = "cache";
	std::string ext = ".ply";
	std::string isofile = bakedir + "/" + framestr + ext;
	exportMeshToFile(isoMesh, isofile);
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

	std::vector<Point3i> extrapolationIndices;
	extrapolationIndices.reserve((isize+1) * jsize * ksize);
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

	std::vector<Point3i> extrapolationIndices;
	extrapolationIndices.reserve(isize * (jsize+1) * ksize);
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

	std::vector<Point3i> extrapolationIndices;
	extrapolationIndices.reserve(isize * jsize * (ksize + 1));
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

void FluidSimulation::advectVelocityField() {
	advectVelocityFieldU();
	advectVelocityFieldV();
	advectVelocityFieldW();
}

void FluidSimulation::applyConstantBodyForces(double dt) {
	Vector3f totalBodyForce = getConstantBodyForce();

	if (fabs(totalBodyForce.x) > 0.0) {
		for (int k = 0; k < ksize; k++) {
			for (int j = 0; j < jsize; j++) {
				for (int i = 0; i < isize + 1; i++) {
					if (materialGrid.isFaceBorderingMaterialU(i, j, k, Material::fluid)) {
						macGrid.setU(i, j, k, macGrid.U(i, j, k) + totalBodyForce.x * dt);
					}
				}
			}
		}

		for (int k = 0; k < ksize; k++) {
			for (int j = 0; j < jsize + 1; j++) {
				for (int i = 0; i < isize; i++) {
					if (materialGrid.isFaceBorderingMaterialV(i, j, k, Material::fluid)) {
						macGrid.setV(i, j, k, macGrid.V(i, j, k) + totalBodyForce.y * dt);
					}
				}
			}
		}

		for (int k = 0; k < ksize + 1; k++) {
			for (int j = 0; j < jsize; j++) {
				for (int i = 0; i < isize; i++) {
					if (materialGrid.isFaceBorderingMaterialW(i, j, k, Material::fluid)) {
						macGrid.setW(i, j, k, macGrid.U(i, j, k) + totalBodyForce.z * dt);
					}
				}
			}
		}
	}
}

void FluidSimulation::applyBodyForcesToVelocityField(double dt) {
	applyConstantBodyForces(dt);
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

void FluidSimulation::applyPressureToFaceU(int i, int j, int k, Array3D<float>& pressureGrid, MACGrid& updatedMACGrid, double deltaTime) {
	double usolid = 0.0;
	double scale = deltaTime / (density * dcell);
	double invscale = 1.0 / scale;

	double p0, p1;
	if (!materialGrid.isCellSolid(i - 1, j, k) && !materialGrid.isCellSolid(i, j, k)) {
		p0 = pressureGrid(i - 1, j, k);
		p1 = pressureGrid(i, j, k);
	} else if (materialGrid.isCellSolid(i - 1, j, k)) {
		p0 = pressureGrid(i, j, k) - invscale * (macGrid.U(i, j, k) - usolid);
		p1 = pressureGrid(i, j, k);
	} else {
		p0 = pressureGrid(i - 1, j, k);
		p1 = pressureGrid(i - 1, j, k) - invscale * (macGrid.U(i, j, k) - usolid);
	}

	double u = macGrid.U(i, j, k) - scale * (p1 - p0);
	updatedMACGrid.setU(i, j, k, u);
}

void FluidSimulation::applyPressureToFaceV(int i, int j, int k, Array3D<float>& pressureGrid, MACGrid& updatedMACGrid, double deltaTime) {
	double vsolid = 0.0;
	double scale = deltaTime / (density * dcell);
	double invscale = 1.0 / scale;

	double p0, p1;
	if (!materialGrid.isCellSolid(i, j - 1, k) && !materialGrid.isCellSolid(i, j, k)) {
		p0 = pressureGrid(i, j - 1, k);
		p1 = pressureGrid(i, j, k);
	} else if (materialGrid.isCellSolid(i, j - 1, k)) {
		p0 = pressureGrid(i, j, k) - invscale * (macGrid.V(i, j, k) - vsolid);
		p1 = pressureGrid(i, j, k);
	} else {
		p0 = pressureGrid(i, j - 1, k);
		p1 = pressureGrid(i, j - 1, k) - invscale * (macGrid.V(i, j, k) - vsolid);
	}

	double v = macGrid.V(i, j, k) - scale * (p1 - p0);
	updatedMACGrid.setV(i, j, k, v);
}

void FluidSimulation::applyPressureToFaceW(int i, int j, int k, Array3D<float>& pressureGrid, MACGrid& updatedMACGrid, double deltaTime) {
	double wsolid = 0.0;
	double scale = deltaTime / (density * dcell);
	double invscale = 1.0 / scale;

	double p0, p1;
	if (!materialGrid.isCellSolid(i, j, k - 1) && !materialGrid.isCellSolid(i, j, k)) {
		p0 = pressureGrid(i, j, k - 1);
		p1 = pressureGrid(i, j, k);
	} else if (materialGrid.isCellSolid(i, j, k - 1)) {
		p0 = pressureGrid(i, j, k) - invscale * (macGrid.W(i, j, k) - wsolid);
		p1 = pressureGrid(i, j, k);
	} else {
		p0 = pressureGrid(i, j, k - 1);
		p1 = pressureGrid(i, j, k - 1) - invscale * (macGrid.W(i, j, k) - wsolid);
	}

	double w = macGrid.W(i, j, k) - scale * (p1 - p0);
	updatedMACGrid.setW(i, j, k, w);
}

void FluidSimulation::applyPressureToVelocityField(Array3D<float>& pressureGrid, double dt) {
	MACGrid updatedMACGrid(isize, jsize, ksize, dcell);

	// TODO: potential for optimization
	// calculate du, dv, dw values
	for (int k = 0; k < ksize; k++) {
		for (int j = 0; j < jsize; j++) {
			for (int i = 0; i < isize + 1; i++) {
				if (materialGrid.isFaceBorderingMaterialU(i, j, k, Material::solid)) {
					updatedMACGrid.setU(i, j, k, 0.0);
				}

				if (materialGrid.isFaceBorderingMaterialU(i, j, k, Material::fluid) && materialGrid.isFaceBorderingMaterialU(i, j, k, Material::solid)) {
					applyPressureToFaceU(i, j, k, pressureGrid, updatedMACGrid, dt);
				}
			}
		}
	}

	for (int k = 0; k < ksize; k++) {
		for (int j = 0; j < jsize + 1; j++) {
			for (int i = 0; i < isize; i++) {
				if (materialGrid.isFaceBorderingMaterialV(i, j, k, Material::solid)) {
					updatedMACGrid.setV(i, j, k, 0.0);
				}

				if (materialGrid.isFaceBorderingMaterialV(i, j, k, Material::fluid) && materialGrid.isFaceBorderingMaterialV(i, j, k, Material::solid)) {
					applyPressureToFaceV(i, j, k, pressureGrid, updatedMACGrid, dt);
				}
			}
		}
	}

	for (int k = 0; k < ksize + 1; k++) {
		for (int j = 0; j < jsize; j++) {
			for (int i = 0; i < isize; i++) {
				if (materialGrid.isFaceBorderingMaterialW(i, j, k, Material::solid)) {
					updatedMACGrid.setW(i, j, k, 0.0);
				}

				if (materialGrid.isFaceBorderingMaterialW(i, j, k, Material::fluid) && materialGrid.isFaceBorderingMaterialW(i, j, k, Material::solid)) {
					applyPressureToFaceW(i, j, k, pressureGrid, updatedMACGrid, dt);
				}
			}
		}
	}

	// add du, dv, dw values back to grid
	for (int k = 0; k < ksize; k++) {
		for (int j = 0; j < jsize; j++) {
			for (int i = 0; i < isize + 1; i++) {
				if (materialGrid.isFaceBorderingMaterialU(i, j, k, Material::fluid)) {
					macGrid.setU(i, j, k, updatedMACGrid.U(i, j, k));
				}
			}
		}
	}

	for (int k = 0; k < ksize; k++) {
		for (int j = 0; j < jsize + 1; j++) {
			for (int i = 0; i < isize; i++) {
				if (materialGrid.isFaceBorderingMaterialV(i, j, k, Material::fluid)) {
					macGrid.setV(i, j, k, updatedMACGrid.V(i, j, k));
				}
			}
		}
	}

	for (int k = 0; k < ksize + 1; k++) {
		for (int j = 0; j < jsize; j++) {
			for (int i = 0; i < isize; i++) {
				if (materialGrid.isFaceBorderingMaterialW(i, j, k, Material::fluid)) {
					macGrid.setW(i, j, k, updatedMACGrid.W(i, j, k));
				}
			}
		}
	}
}

void FluidSimulation::extrapolateFluidVelocities(MACGrid& velocityGrid) {
	int numLayers = (int)std::ceil(CFLConditionNumber + 2);
	velocityGrid.extrapolateVelocityField(materialGrid, numLayers);
}

void FluidSimulation::updateMarkerParticleVelocitiesSubset(int start, int end) {
	int size = end - start + 1;
	std::vector<Point3f> positions;
	positions.reserve(size);
	std::vector<Vector3f> velocityNew, velocityOld;
	velocityNew.reserve(size);
	velocityOld.reserve(size);

	for (int idx = start; idx <= end; idx++) {
		positions.push_back(markerParticles[idx].position);
	}

	particleAdvector.tricubicInterpolate(positions, &macGrid, velocityNew);
	particleAdvector.tricubicInterpolate(positions, &prevMACGrid, velocityOld);

	for (int i = 0; i < positions.size(); i++) {
		MarkerParticle mp = markerParticles[start + i];

		Vector3f vPIC = velocityNew[i];
		Vector3f vFLIP = mp.velocity + velocityNew[i] - velocityOld[i];

		Vector3f v = vFLIP * ratioFLIPPIC + vPIC * (1 - ratioFLIPPIC);
		markerParticles[start + i].velocity = v;
	}
}

void FluidSimulation::updateMarkerParticleVelocities() {
	for (int start = 0; start < markerParticles.size(); start += maxParticlesPerFLIPPICUpdate) {
		int end = start + maxParticlesPerFLIPPICUpdate - 1;
		end = fmin(end, markerParticles.size() - 1);

		updateMarkerParticleVelocitiesSubset(start, end);
	}
}

void FluidSimulation::advanceMarkerParticlesSubset(int start, int end, double dt) {
	std::vector<Point3f> positions;
	positions.reserve(end - start + 1);
	for (int i = start; i <= end; i++) {
		positions.push_back(markerParticles[i].position);
	}

	std::vector<Point3f> output;
	particleAdvector.advectParticles(positions, &macGrid, dt, output);

	for (int i = 0; i < output.size(); i++) {
		MarkerParticle mp = markerParticles[start + i];
		Point3f nextpos = output[i];

		Point3i idx = positionToGridIndex(nextpos.x, nextpos.y, nextpos.z, dcell);
		if (materialGrid.isCellSolid(idx.x, idx.y, idx.z)) {
			nextpos = resolveParticleSolidCollision(mp.position, nextpos);
		}

		markerParticles[start + i].position = nextpos;
	}
}

void FluidSimulation::removeMarkerParticles() {
	for (int i = markerParticles.size() - 2; i >= 0; i--) {
		int j = (rand() % (i - 0 + 1));
		MarkerParticle tmp = markerParticles[i];
		markerParticles[i] = markerParticles[j];
		markerParticles[j] = tmp;
	}

	std::vector<bool> isRemoved(markerParticles.size());
	Array3D<int> countGrid(isize, jsize, ksize, 0);
	for (int i = 0; i < markerParticles.size(); i++) {
		MarkerParticle mp = markerParticles[i];
		Point3i idx = positionToGridIndex(mp.position.x, mp.position.y, mp.position.z, dcell);
		if (countGrid(idx.x, idx.y, idx.z) >= maxMarkerParticlesPerCell) {
			isRemoved.push_back(true);
			continue;
		}
		countGrid.set(idx.x, idx.y, idx.z, countGrid(idx.x, idx.y, idx.z) + 1);
		isRemoved.push_back(false);
	}

	markerParticles.erase(std::remove_if(markerParticles.begin(), markerParticles.end(),
								  [&isRemoved, &markerParticles = markerParticles](auto const& i) {
									  return isRemoved.at(&i - markerParticles.data());
								  }),
			markerParticles.end());
	markerParticles.shrink_to_fit();
}

Point3f FluidSimulation::resolveParticleSolidCollision(Point3f oldpos, Point3f newpos) {
	Point3i oldidx = positionToGridIndex(oldpos.x, oldpos.y, oldpos.z, dcell);
	Point3i newidx = positionToGridIndex(newpos.x, newpos.y, newpos.z, dcell);
	if (!materialGrid.isCellSolid(oldidx.x, oldidx.y, oldidx.z) || materialGrid.isCellSolid(newidx.x, newidx.y, newidx.z)) {
		return oldpos;
	}

	Point3i voxel;
	bool foundVoxel = getLineSegmentVoxelIntersection(oldpos, newpos, dcell, materialGrid, &voxel);

	if (!foundVoxel) {
		return oldpos;
	}

	Vector3f line = normalize(newpos - oldpos);
	Point3f voxelpos = gridIndexToPosition(voxel.x, voxel.y, voxel.z, dcell);
	Bounds3f bbox(voxelpos, Point3f(voxelpos.x + dcell, voxelpos.y + dcell, voxelpos.z + dcell));

	float t0, t1;
	if (!bbox.intersectP(oldpos, line, Infinity, &t0, &t1)) {
		return oldpos;
	}

	Point3f isectP = oldpos + line * t0;
	Point3f resolvedP = isectP - line * 0.05f * dcell;

	Point3i resolvedidx = positionToGridIndex(resolvedP.x, resolvedP.y, resolvedP.z, dcell);
	if (materialGrid.isCellSolid(resolvedidx.x, resolvedidx.y, resolvedidx.z)) {
		return oldpos;
	}

	return resolvedP;
}

void FluidSimulation::advanceMarkerParticles(double dt) {
	for (int start = 0; start < markerParticles.size(); start += maxParticlesPerFLIPPICUpdate) {
		int end = start + maxParticlesPerFLIPPICUpdate - 1;
		end = fmin(end, markerParticles.size() - 1);

		advanceMarkerParticlesSubset(start, end, dt);
	}

	removeMarkerParticles();
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
	reconstructFluidSurfaceMesh();

	// advect velocity field with new marker particles
	advectVelocityField();
	prevMACGrid = macGrid;
	extrapolateFluidVelocities(prevMACGrid);

	// apply body forces, e.g. gravity
	applyBodyForcesToVelocityField(dt);

	// update and apply pressure
	Array3D<float> pressureGrid(isize, jsize, ksize, 0.0f);
	updatePressureGrid(pressureGrid, dt);
	applyPressureToVelocityField(pressureGrid, dt);

	extrapolateFluidVelocities(macGrid);

	// update and advect marker particles
	updateMarkerParticleVelocities();
	advanceMarkerParticles(dt);
}

} // namespace foc