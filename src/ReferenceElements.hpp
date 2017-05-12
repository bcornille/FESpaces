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

#endif
