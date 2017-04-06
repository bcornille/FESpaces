#include "Basis1D.hpp"

#ifndef _ReferenceElements_hpp
#define _ReferenceElements_hpp

using namespace Eigen;

class H1_1D
{
	public:
		H1_1D(int p = 0);
		~H1_1D() = default;
		VectorXd eval(double x);
		VectorXd evalD(double x);
		int dofs();
	private:
		GaussLobatto gll;
};

H1_1D::H1_1D(int p) : gll(p + 1) {}

inline VectorXd H1_1D::eval(double x)
{
	return gll.evalGLL(x).leftCols<1>();
}

inline VectorXd H1_1D::evalD(double x)
{
	return gll.evalGLL(x).rightCols<1>();
}

inline int H1_1D::dofs()
{
	return gll.getN();
}

class L2_1D
{
	public:
		L2_1D(in p = 0);
		~L2_1D() = default;
		VectorXd eval(double x);
		int dofs();
	private:
		GaussLegendre gl;
};

L2_1D::L2_1D(int p) : gl(p) {}

inline VectorXd L2_1D::eval(double x)
{
	return gl.evalGL(x).leftCols<1>();
}

inline int L2_1D::dofs()
{
	return gl.getN();
}

#endif
