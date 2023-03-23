#include "pressuresolver.h"

namespace foc {

VectorXd::VectorXd() {
}

VectorXd::VectorXd(int size) :
		vector(size, 0.0) {
}

VectorXd::VectorXd(int size, double fill) :
		vector(size, fill) {
}

VectorXd::VectorXd(VectorXd& v) {
	vector.reserve(v.size());
	for (unsigned int i = 0; i < v.size(); i++) {
		vector.push_back(v[i]);
	}
}

VectorXd::~VectorXd() {
}

const double VectorXd::operator[](int i) const {
	return vector[i];
}

double& VectorXd::operator[](int i) {
	return vector[i];
}

void VectorXd::fill(double fill) {
	for (unsigned int i = 0; i < vector.size(); i++) {
		vector[i] = fill;
	}
}

double VectorXd::dot(VectorXd& v) {
	double sum = 0.0;
	for (unsigned int i = 0; i < vector.size(); i++) {
		sum += vector[i] * v.vector[i];
	}

	return sum;
}

double VectorXd::absMaxCoeff() {
	double max = -std::numeric_limits<double>::infinity();
	for (unsigned int i = 0; i < vector.size(); i++) {
		if (fabs(vector[i]) > max) {
			max = fabs(vector[i]);
		}
	}

	return max;
}

MatrixCoefficients::MatrixCoefficients() {
}

MatrixCoefficients::MatrixCoefficients(int size) :
		cells(size, MatrixCell()) {
}

MatrixCoefficients::~MatrixCoefficients() {
}

const MatrixCell MatrixCoefficients::operator[](int i) const {
	return cells[i];
}

MatrixCell& MatrixCoefficients::operator[](int i) {
	return cells[i];
}

PressureSolver::PressureSolver() {
}

PressureSolver::~PressureSolver() {
}

void PressureSolver::solve(PressureSolverParameters params, VectorXd& pressure) {
	initialize(params);

	pressure.fill(0.0);

	keymap.reserve(isize * jsize * ksize);
	for (int i = 0; i < fluidCells->size(); i++) {
		Point3i idx = fluidCells->at(i);
		keymap[getFlatIndex(idx.x, idx.y, idx.z)] = i;
	}

	VectorXd b(matSize);
	calculateNegativeDivergenceVector(b);
	if (b.absMaxCoeff() < pressureSolveTolerance) {
		return;
	}

	MatrixCoefficients A(matSize);
	calculateMatrixCoefficients(A);

	VectorXd precon(matSize);
	calculatePreconditionerVector(A, precon);

	solvePressureSystem(A, b, precon, pressure);
}

void PressureSolver::initialize(PressureSolverParameters params) {
	isize = params.materialGrid->isize;
	jsize = params.materialGrid->jsize;
	ksize = params.materialGrid->ksize;
	cellsize = params.cellsize;
	density = params.density;
	deltaTime = params.deltaTime;
	fluidCells = params.fluidCells;
	materialGrid = params.materialGrid;
	vField = params.velocityField;
	matSize = (int)fluidCells->size();
}

void PressureSolver::calculateNegativeDivergenceVector(VectorXd& b) {
	double scale = 1.0 / cellsize;
	for (int idx = 0; idx < fluidCells->size(); idx++) {
		int i = fluidCells->at(idx).x;
		int j = fluidCells->at(idx).y;
		int k = fluidCells->at(idx).z;

		double value = -scale * (double)(vField->U(i + 1, j, k) - vField->U(i, j, k) + vField->V(i, j + 1, k) - vField->V(i, j, k) + vField->W(i, j, k + 1) - vField->W(i, j, k));

		int vecidx = getFlatIndex(i, j, k);
		b[keymap[vecidx]] = value;
	}

	float usolid = 0.0;
	float vsolid = 0.0;
	float wsolid = 0.0;
	for (int idx = 0; idx < fluidCells->size(); idx++) {
		int i = fluidCells->at(idx).x;
		int j = fluidCells->at(idx).y;
		int k = fluidCells->at(idx).z;
		int vecidx = keymap[getFlatIndex(i, j, k)];

		if (materialGrid->isCellSolid(i - 1, j, k)) {
			b[vecidx] -= (float)scale * (vField->U(i, j, k) - usolid);
		}
		if (materialGrid->isCellSolid(i + 1, j, k)) {
			b[vecidx] += (float)scale * (vField->U(i + 1, j, k) - usolid);
		}

		if (materialGrid->isCellSolid(i, j - 1, k)) {
			b[vecidx] -= (float)scale * (vField->V(i, j, k) - vsolid);
		}
		if (materialGrid->isCellSolid(i, j + 1, k)) {
			b[vecidx] += (float)scale * (vField->V(i, j + 1, k) - vsolid);
		}

		if (materialGrid->isCellSolid(i, j, k - 1)) {
			b[vecidx] -= (float)scale * (vField->W(i, j, k) - wsolid);
		}
		if (materialGrid->isCellSolid(i, j, k + 1)) {
			b[vecidx] += (float)scale * (vField->W(i, j, k + 1) - wsolid);
		}
	}
}

int PressureSolver::getNumFluidOrAirCellNeighbours(int i, int j, int k) {
	int n = 0;
	if (!materialGrid->isCellSolid(i - 1, j, k)) {
		n++;
	}
	if (!materialGrid->isCellSolid(i + 1, j, k)) {
		n++;
	}
	if (!materialGrid->isCellSolid(i, j - 1, k)) {
		n++;
	}
	if (!materialGrid->isCellSolid(i, j + 1, k)) {
		n++;
	}
	if (!materialGrid->isCellSolid(i, j, k - 1)) {
		n++;
	}
	if (!materialGrid->isCellSolid(i, j, k + 1)) {
		n++;
	}

	return n;
}

void PressureSolver::calculateMatrixCoefficients(MatrixCoefficients& A) {
	for (unsigned int idx = 0; idx < fluidCells->size(); idx++) {
		int i = fluidCells->at(idx).x;
		int j = fluidCells->at(idx).y;
		int k = fluidCells->at(idx).z;
		int vecidx = keymap[getFlatIndex(i, j, k)];

		int n = getNumFluidOrAirCellNeighbours(i, j, k);
		A.cells[vecidx].diag = (char)n;

		if (materialGrid->isCellFluid(i + 1, j, k)) {
			A.cells[vecidx].plusi = 0x01;
		}

		if (materialGrid->isCellFluid(i, j + 1, k)) {
			A.cells[vecidx].plusj = 0x01;
		}

		if (materialGrid->isCellFluid(i, j, k + 1)) {
			A.cells[vecidx].plusk = 0x01;
		}
	}
}

void PressureSolver::calculatePreconditionerVector(MatrixCoefficients& A, VectorXd& precon) {
	double scale = deltaTime / (density * cellsize * cellsize);
	double negscale = -scale;

	double tau = 0.97; // Tuning constant
	double sigma = 0.25; // safety constant
	for (unsigned int idx = 0; idx < fluidCells->size(); idx++) {
		int i = fluidCells->at(idx).x;
		int j = fluidCells->at(idx).y;
		int k = fluidCells->at(idx).z;
		int vecidx = keymap[getFlatIndex(i, j, k)];

		int vecidx_im1 = keymap[getFlatIndex(i - 1, j, k)];
		int vecidx_jm1 = keymap[getFlatIndex(i, j - 1, k)];
		int vecidx_km1 = keymap[getFlatIndex(i, j, k - 1)];

		double diag = (double)A[vecidx].diag * scale;

		double plusi_im1 = vecidx_im1 != -1 ? (double)A[vecidx_im1].plusi * negscale : 0.0;
		double plusi_jm1 = vecidx_jm1 != -1 ? (double)A[vecidx_jm1].plusi * negscale : 0.0;
		double plusi_km1 = vecidx_km1 != -1 ? (double)A[vecidx_km1].plusi * negscale : 0.0;

		double plusj_im1 = vecidx_im1 != -1 ? (double)A[vecidx_im1].plusj * negscale : 0.0;
		double plusj_jm1 = vecidx_jm1 != -1 ? (double)A[vecidx_jm1].plusj * negscale : 0.0;
		double plusj_km1 = vecidx_km1 != -1 ? (double)A[vecidx_km1].plusj * negscale : 0.0;

		double plusk_im1 = vecidx_im1 != -1 ? (double)A[vecidx_im1].plusk * negscale : 0.0;
		double plusk_jm1 = vecidx_jm1 != -1 ? (double)A[vecidx_jm1].plusk * negscale : 0.0;
		double plusk_km1 = vecidx_km1 != -1 ? (double)A[vecidx_km1].plusk * negscale : 0.0;

		double precon_im1 = vecidx_im1 != -1 ? precon[vecidx_im1] : 0.0;
		double precon_jm1 = vecidx_jm1 != -1 ? precon[vecidx_jm1] : 0.0;
		double precon_km1 = vecidx_km1 != -1 ? precon[vecidx_km1] : 0.0;

		double v1 = plusi_im1 * precon_im1;
		double v2 = plusj_jm1 * precon_jm1;
		double v3 = plusk_km1 * precon_km1;
		double v4 = precon_im1 * precon_im1;
		double v5 = precon_jm1 * precon_jm1;
		double v6 = precon_km1 * precon_km1;

		double e = diag - v1 * v1 - v2 * v2 - v3 * v3 -
				tau * (plusi_im1 * (plusj_im1 + plusk_im1) * v4 + plusj_jm1 * (plusi_jm1 + plusk_jm1) * v5 + plusk_km1 * (plusi_km1 + plusj_km1) * v6);

		if (e < sigma * diag) {
			e = diag;
		}

		if (fabs(e) > 10e-9) {
			precon[vecidx] = 1.0 / sqrt(e);
		}
	}
}

void PressureSolver::applyPreconditioner(MatrixCoefficients& A, VectorXd& precon, VectorXd& residual, VectorXd& vect) {
	double scale = deltaTime / (density * cellsize * cellsize);
	double negscale = -scale;

	// Solve A*q = residual
	VectorXd q(matSize);
	for (unsigned int idx = 0; idx < fluidCells->size(); idx++) {
		int i = fluidCells->at(idx).x;
		int j = fluidCells->at(idx).y;
		int k = fluidCells->at(idx).z;
		int vecidx = keymap[getFlatIndex(i, j, k)];

		int vecidx_im1 = keymap[getFlatIndex(i - 1, j, k)];
		int vecidx_jm1 = keymap[getFlatIndex(i, j - 1, k)];
		int vecidx_km1 = keymap[getFlatIndex(i, j, k - 1)];

		double plusi_im1 = 0.0;
		double precon_im1 = 0.0;
		double q_im1 = 0.0;
		if (vecidx_im1 != -1) {
			plusi_im1  = (double)A[vecidx_im1].plusi * negscale;
			precon_im1 = precon[vecidx_im1];
			q_im1      = q[vecidx_im1];
		}

		double plusj_jm1 = 0.0;
		double precon_jm1 = 0.0;
		double q_jm1 = 0.0;
		if (vecidx_jm1 != -1) {
			plusj_jm1  = (double)A[vecidx_jm1].plusj * negscale;
			precon_jm1 = precon[vecidx_jm1];
			q_jm1      = q[vecidx_jm1];
		}

		double plusk_km1 = 0.0;
		double precon_km1 = 0.0;
		double q_km1 = 0.0;
		if (vecidx_km1 != -1) {
			plusk_km1  = (double)A[vecidx_km1].plusk * negscale;
			precon_km1 = precon[vecidx_km1];
			q_km1      = q[vecidx_km1];
		}

		double t = residual[vecidx] - plusi_im1 * precon_im1 * q_im1 -
				plusj_jm1 * precon_jm1 * q_jm1 -
				plusk_km1 * precon_km1 * q_km1;

		t = t*precon[vecidx];
		q[vecidx] = t;
	}

	// Solve transpose(A)*z = q
	for (int idx = (int)fluidCells->size() - 1; idx >= 0; idx--) {
		int i = fluidCells->at(idx).x;
		int j = fluidCells->at(idx).y;
		int k = fluidCells->at(idx).z;
		int vecidx = keymap[getFlatIndex(i, j, k)];

		int vecidx_ip1 = keymap[getFlatIndex(i + 1, j, k)];
		int vecidx_jp1 = keymap[getFlatIndex(i, j + 1, k)];
		int vecidx_kp1 = keymap[getFlatIndex(i, j, k + 1)];

		double vect_ip1 = vecidx_ip1 != -1 ? vect[vecidx_ip1] : 0.0;
		double vect_jp1 = vecidx_jp1 != -1 ? vect[vecidx_jp1] : 0.0;
		double vect_kp1 = vecidx_kp1 != -1 ? vect[vecidx_kp1] : 0.0;

		double plusi = (double)A[vecidx].plusi * negscale;
		double plusj = (double)A[vecidx].plusj * negscale;
		double plusk = (double)A[vecidx].plusk * negscale;

		double preconval = precon[vecidx];
		double t = q[vecidx] - plusi * preconval * vect_ip1 -
				plusj * preconval * vect_jp1 -
				plusk * preconval * vect_kp1;

		t = t*preconval;
		vect[vecidx] = t;
	}
}

void PressureSolver::applyMatrix(MatrixCoefficients& A, VectorXd& x, VectorXd& result) {
	double scale = deltaTime / (density * cellsize * cellsize);
	double negscale = -scale;

	for (unsigned int idx = 0; idx < fluidCells->size(); idx++) {
		int i = fluidCells->at(idx).x;
		int j = fluidCells->at(idx).y;
		int k = fluidCells->at(idx).z;

		// val = dot product of column vector x and idxth row of matrix A
		double val = 0.0;
		int vecidx = keymap[getFlatIndex(i - 1, j, k)];
		if (vecidx != -1) { val += x.vector[vecidx]; }

		vecidx = keymap[getFlatIndex(i + 1, j, k)];
		if (vecidx != -1) { val += x.vector[vecidx]; }

		vecidx = keymap[getFlatIndex(i, j - 1, k)];
		if (vecidx != -1) { val += x.vector[vecidx]; }

		vecidx = keymap[getFlatIndex(i, j + 1, k)];
		if (vecidx != -1) { val += x.vector[vecidx]; }

		vecidx = keymap[getFlatIndex(i, j, k - 1)];
		if (vecidx != -1) { val += x.vector[vecidx]; }

		vecidx = keymap[getFlatIndex(i, j, k + 1)];
		if (vecidx != -1) { val += x.vector[vecidx]; }

		val *= negscale;

		vecidx = keymap[getFlatIndex(i, j, k)];
		val += (double)A.cells[vecidx].diag * scale * x.vector[vecidx];

		result.vector[vecidx] = val;
	}
}

// v1 += v2*scale
void PressureSolver::addScaledVector(VectorXd& v1, VectorXd& v2, double scale) {
	for (unsigned int idx = 0; idx < v1.size(); idx++) {
		v1.vector[idx] += v2.vector[idx]*scale;
	}
}

// result = v1*s1 + v2*s2
void PressureSolver::addScaledVectors(VectorXd& v1, double s1, VectorXd& v2, double s2, VectorXd& result) {
	for (unsigned int idx = 0; idx < v1.size(); idx++) {
		result.vector[idx] = v1.vector[idx]*s1 + v2.vector[idx]*s2;
	}
}

// Solve (A*pressure = b) with Modified Incomplete Cholesky
// Conjugate Gradient method (MICCG(0))
void PressureSolver::solvePressureSystem(MatrixCoefficients& A, VectorXd& b, VectorXd& precon, VectorXd& pressure) {
	double tol = pressureSolveTolerance;
	if (b.absMaxCoeff() < tol) {
		return;
	}

	VectorXd residual(b);
	VectorXd auxillary(matSize);
	applyPreconditioner(A, precon, residual, auxillary);

	VectorXd search(auxillary);

	double alpha = 0.0;
	double beta = 0.0;
	double sigma = auxillary.dot(residual);
	double sigmaNew = 0.0;
	int iterationNumber = 0;

	while (iterationNumber < maxCGIterations) {
		applyMatrix(A, search, auxillary);
		alpha = sigma / auxillary.dot(search);
		addScaledVector(pressure, search, alpha);
		addScaledVector(residual, auxillary, -alpha);

		if (residual.absMaxCoeff() < tol) {
			return;
		}

		applyPreconditioner(A, precon, residual, auxillary);
		sigmaNew = auxillary.dot(residual);
		beta = sigmaNew / sigma;
		addScaledVectors(auxillary, 1.0, search, beta, search);
		sigma = sigmaNew;

		iterationNumber++;
	}
}

} // namespace foc