#include "Basis1D.hpp"

#ifndef _ReferenceElements_hpp
#define _ReferenceElements_hpp

using namespace Eigen;

class H1_1D
{
	public:
		H1_1D(int p = 1);
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
		L2_1D(int p = 1);
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

class L2_1D_EF
{
	public:
		L2_1D_EF(int p = 1);
		~L2_1D_EF() = default;
		VectorXd eval(double x);
		int dofs();
	private:
		EdgeFunction ef;
};

L2_1D_EF::L2_1D_EF(int p) : ef(p) {}

inline VectorXd L2_1D_EF::eval(double x)
{
	return ef.evalEF(x);
}

inline int L2_1D_EF::dofs()
{
	return ef.getN();
}

class H1_2D
{
	public:
		H1_2D(int p = 1);
		~H1_2D() = default;
		VectorXd eval(Vector2d x);
		MatrixX2d evalGrad(Vector2d x);
		int dofs();
	private:
		int n_dof;
		GaussLobatto gll;
};

H1_2D::H1_2D(int p) : n_dof((p + 1)*(p + 1)), gll(p + 1) {}

inline VectorXd H1_2D::eval(Vector2d x)
{
	VectorXd vals(n_dof);
	MatrixXd mat_vals = (gll.evalGLL(x[0]).leftCols<1>()
		*gll.evalGLL(x[1]).leftCols<1>().transpose());
	vals = Map<VectorXd>(mat_vals.data(), n_dof);
	return vals;
}

inline MatrixX2d H1_2D::evalGrad(Vector2d x)
{
	MatrixX2d grad(n_dof, 2);
	MatrixX2d gll_x = gll.evalGLL(x[0]);
	MatrixX2d gll_y = gll.evalGLL(x[1]);
	MatrixXd mat_vals = gll_x.rightCols<1>()*gll_y.leftCols<1>().transpose();
	grad.leftCols<1>() = Map<VectorXd>(mat_vals.data(), n_dof);
	mat_vals = gll_x.leftCols<1>()*gll_y.rightCols<1>().transpose();
	grad.rightCols<1>() = Map<VectorXd>(mat_vals.data(), n_dof);
	return grad;
}

inline int H1_2D::dofs()
{
	return n_dof;
}

class HCurl_2D
{
	public:
		HCurl_2D(int p = 1);
		~HCurl_2D() = default;
		MatrixX2d eval(Vector2d x);
		int dofs();
	private:
		int n_dof;
		GaussLobatto gll;
		GaussLegendre gl;
};

HCurl_2D::HCurl_2D(int p) : n_dof(p*(p + 1)), gll(p + 1), gl(p) {}

inline MatrixX2d HCurl_2D::eval(Vector2d x)
{
	MatrixX2d vals(n_dof, 2);
	MatrixXd mat_xvals = (gl.evalGL(x[0]).leftCols<1>()
		*gll.evalGLL(x[1]).leftCols<1>().transpose());
	vals.leftCols<1>() = Map<VectorXd>(mat_xvals.data(), n_dof);
	MatrixXd mat_yvals = (gll.evalGLL(x[0]).leftCols<1>()
		*gl.evalGL(x[1]).leftCols<1>().transpose());
	vals.rightCols<1>() = Map<VectorXd>(mat_yvals.data(), n_dof);
	return vals;
}

inline int HCurl_2D::dofs()
{
	return n_dof;
}

#endif
