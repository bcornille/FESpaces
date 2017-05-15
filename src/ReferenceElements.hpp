#include "Basis1D.hpp"

#ifndef _ReferenceElements_hpp
#define _ReferenceElements_hpp

using namespace Eigen;

/**
 * @brief      Class for $H^1$ elements on $[-1, 1]$.
 */
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

/**
 * @brief      Constructs the object.
 *
 * @param[in]  p     Polynomial degree.
 */
H1_1D::H1_1D(int p) : gll(p + 1) {}

/**
 * @brief      Evaluates the basis functions at the desired point.
 *
 * @param[in]  x     $x \in [-1, 1]$.
 *
 * @return     Vector of basis function evaluations.
 */
inline VectorXd H1_1D::eval(double x)
{
	return gll.evalGLL(x).leftCols<1>();
}

/**
 * @brief      Evaluates the derivatives of the basis functions at the desired
 *             point.
 *
 * @param[in]  x     $x \in [-1, 1]$.
 *
 * @return     Vector of the basis function derivative evaluations.
 */
inline VectorXd H1_1D::evalD(double x)
{
	return gll.evalGLL(x).rightCols<1>();
}

/**
 * @brief      Get number of degrees of freedom.
 *
 * @return     Number of DoF.
 */
inline int H1_1D::dofs()
{
	return gll.getN();
}

/**
 * @brief      Class for $L^2$ element on $[-1, 1]$.
 */
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

/**
 * @brief      Constructs the object.
 *
 * @param[in]  p     Number of degrees of freedom, or polynomial degree plus 1.
 */
L2_1D::L2_1D(int p) : gl(p) {}

/**
 * @brief      Evaluate the basis functions at the desired point.
 *
 * @param[in]  x     $x \in [-1, 1]$.
 *
 * @return     Vector of the basis function evaluations.
 */
inline VectorXd L2_1D::eval(double x)
{
	return gl.evalGL(x).leftCols<1>();
}

/**
 * @brief      Get number of degrees of freedom.
 *
 * @return     Number of DoF.
 */
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
		int dof_per_dim;
		int n_dof;
		GaussLobatto gll;
		GaussLegendre gl;
};

HCurl_2D::HCurl_2D(int p) : dof_per_dim(p*(p + 1)), n_dof(2*p*(p + 1)), gll(p + 1), gl(p) {}

inline MatrixX2d HCurl_2D::eval(Vector2d x)
{
	MatrixX2d vals(n_dof, 2);
	vals.setZero();
	MatrixXd mat_xvals = (gl.evalGL(x[0]).leftCols<1>()
		*gll.evalGLL(x[1]).leftCols<1>().transpose());
	vals.topLeftCorner(dof_per_dim, 1) = Map<VectorXd>(mat_xvals.data(), dof_per_dim);
	MatrixXd mat_yvals = (gll.evalGLL(x[0]).leftCols<1>()
		*gl.evalGL(x[1]).leftCols<1>().transpose());
	vals.bottomRightCorner(dof_per_dim,1) = Map<VectorXd>(mat_yvals.data(), dof_per_dim);
	return vals;
}

inline int HCurl_2D::dofs()
{
	return n_dof;
}

class HDiv_2D
{
	public:
		HDiv_2D(int p = 1);
		~HDiv_2D() = default;
		MatrixX2d eval(Vector2d x);
		VectorXd evalDiv(Vector2d x);
		int dofs();
	private:
		int dof_per_dim;
		int n_dof;
		GaussLobatto gll;
		GaussLegendre gl;
};

HDiv_2D::HDiv_2D(int p) : dof_per_dim((p + 1)*p), n_dof(2*(p + 1)*p), gll(p + 1), gl(p) {}

inline MatrixX2d HDiv_2D::eval(Vector2d x)
{
	MatrixX2d vals(n_dof, 2);
	vals.setZero();
	MatrixXd mat_yvals = (gl.evalGL(x[0]).leftCols<1>()
		*gll.evalGLL(x[1]).leftCols<1>().transpose());
	vals.topRightCorner(dof_per_dim, 1) = Map<VectorXd>(mat_yvals.data(), dof_per_dim);
	MatrixXd mat_xvals = (gll.evalGLL(x[0]).leftCols<1>()
		*gl.evalGL(x[1]).leftCols<1>().transpose());
	vals.bottomLeftCorner(dof_per_dim,1) = Map<VectorXd>(mat_xvals.data(), dof_per_dim);
	return vals;
}

inline VectorXd HDiv_2D::evalDiv(Vector2d x)
{
	VectorXd divg(n_dof);
	MatrixXd mat_yvals = (gl.evalGL(x[0]).leftCols<1>()
		*gll.evalGLL(x[1]).rightCols<1>().transpose());
	divg.head(dof_per_dim) = Map<VectorXd>(mat_yvals.data(), dof_per_dim);
	MatrixXd mat_xvals = (gll.evalGLL(x[0]).rightCols<1>()
		*gl.evalGL(x[1]).leftCols<1>().transpose());
	divg.tail(dof_per_dim) = Map<VectorXd>(mat_xvals.data(), dof_per_dim);
	return divg;
}

inline int HDiv_2D::dofs()
{
	return n_dof;
}

class L2_2D
{
	public:
		L2_2D(int p = 1);
		~L2_2D() = default;
		VectorXd eval(Vector2d x);
		int dofs();
	private:
		int n_dof;
		GaussLegendre gl;
};

L2_2D::L2_2D(int p) : n_dof(p*p), gl(p) {}

inline VectorXd L2_2D::eval(Vector2d x)
{
	VectorXd vals(n_dof);
	MatrixXd mat_vals = (gl.evalGL(x[0]).leftCols<1>()
		*gl.evalGL(x[1]).leftCols<1>().transpose());
	vals = Map<VectorXd>(mat_vals.data(), n_dof);
	return vals;
}

inline int L2_2D::dofs()
{
	return n_dof;
}

#endif
