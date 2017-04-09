#include "json.hpp"
#include "Mesh1D.hpp"
#include "Integrators.hpp"
#include "ReferenceElements.hpp"

#ifndef _PoissonSolver1D_hpp
#define _PoissonSolver1D_hpp

using namespace Eigen;
using namespace nlohmann;

class PoissonSolver1D
{
	public:
		PoissonSolver1D(json params);
		~PoissonSolver1D() = default;
	private:
		SparseMatrix<double> system;
		VectorXd rhs;
		SimplicialLLT<SparseMatrix<double> > solver;
};

PoissonSolver1D::PoissonSolver1D(json params) : 

#endif
