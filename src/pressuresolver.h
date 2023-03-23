#ifndef FOC_PRESSURESOLVER_H
#define FOC_PRESSURESOLVER_H

#include "foc.h"
#include "cellmaterialgrid.h"
#include "macgrid.h"

namespace foc {

struct PressureSolverParameters {
	double cellsize;
	double density;
	double deltaTime;

	std::vector<Point3i>* fluidCells;
	CellMaterialGrid* materialGrid;
	MACGrid* velocityField;
};

class VectorXd {
public:
	VectorXd();
	VectorXd(int size);
	VectorXd(int size, double fill);
	VectorXd(VectorXd& vector);
	~VectorXd();

	const double operator[](int i) const;
	double& operator[](int i);

	inline size_t size() {
		return vector.size();
	}

	void fill(double fill);
	double dot(VectorXd& v);
	double absMaxCoeff();

	std::vector<double> vector;
};

struct MatrixCell {
	char diag;
	char plusi;
	char plusj;
	char plusk;

	MatrixCell() :
			diag(0x00), plusi(0x00), plusj(0x00), plusk(0x00) {}
};

class MatrixCoefficients {
public:
	MatrixCoefficients();
	MatrixCoefficients(int size);
	~MatrixCoefficients();

	const MatrixCell operator[](int i) const;
	MatrixCell& operator[](int i);

	inline size_t size() {
		return cells.size();
	}

	std::vector<MatrixCell> cells;
};

class PressureSolver {
public:
	PressureSolver();
	~PressureSolver();

	void solve(PressureSolverParameters params, VectorXd& pressure);

private:
	void initialize(PressureSolverParameters params);
	void calculateNegativeDivergenceVector(VectorXd& b);
	int getNumFluidOrAirCellNeighbours(int i, int j, int k);
	void calculateMatrixCoefficients(MatrixCoefficients& A);
	void calculatePreconditionerVector(MatrixCoefficients& A, VectorXd& precon);
	void applyPreconditioner(MatrixCoefficients& A,
			VectorXd& precon,
			VectorXd& residual,
			VectorXd& vect);
	void applyMatrix(MatrixCoefficients& A, VectorXd& x, VectorXd& result);
	void addScaledVector(VectorXd& v1, VectorXd& v2, double scale);
	void addScaledVectors(VectorXd& v1, double s1,
			VectorXd& v2, double s2,
			VectorXd& result);
	void solvePressureSystem(MatrixCoefficients& A,
			VectorXd& b,
			VectorXd& precon,
			VectorXd& pressure);


	inline unsigned int getFlatIndex(int i, int j, int k) {
		return (unsigned int)i + (unsigned int)isize *
				((unsigned int)j + (unsigned int)jsize * (unsigned int)k);
	}

private:
	int isize = 0;
	int jsize = 0;
	int ksize = 0;
	double cellsize = 0;
	double density = 0;
	double deltaTime = 0;
	int matSize = 0;

	double pressureSolveTolerance = 1e-6;
	int maxCGIterations = 200;

	std::vector<Point3i>* fluidCells;
	CellMaterialGrid* materialGrid;
	MACGrid* vField;
	std::vector<int> keymap;
};

} // namespace foc

#endif // FOC_PRESSURESOLVER_H
