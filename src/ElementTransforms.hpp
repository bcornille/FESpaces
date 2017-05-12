#include "Eigen/Core"
#include "Basis1D.hpp"

#ifndef _ElementTransforms_hpp
#define _ElementTransforms_hpp

using namespace Eigen;

/**
 * @brief      Class for one-dimensional transforms. Generates a coordinate
 *             transformation map between $[-1,1]$ and a general segment in $x$.
 *             This transformation can be of arbitrary polynomial order.
 */
class Transform1D
{
	public:
		Transform1D(double a = -1.0, double b = 1.0);
		Transform1D(RowVectorXd nodes);
		~Transform1D() = default;
		double forwardTransform(double x_hat);
		double jacobian(double x_hat);
	private:
		RowVectorXd x_nodes;
		GaussLobatto gl;
};

/**
 * @brief      Constructs the object. Creates a linear transformation.
 *
 * @param[in]  a     Minimum $x$ value of the segment.
 * @param[in]  b     Maximum $x$ value of the segment.
 */
Transform1D::Transform1D(double a, double b) : x_nodes(2), gl(2)
{
	x_nodes[0] = a;
	x_nodes[1] = b;
}

/**
 * @brief      Constructs the object. Creates an arbitrary transformation.
 *
 * @param[in]  nodes  The nodes corresponding to the transformed GLL points.
 */
Transform1D::Transform1D(RowVectorXd nodes) : x_nodes(nodes), gl(nodes.size()) {}

/**
 * @brief      Calculates the forward transformation $x=F(\hat{x})$.
 *
 * @param[in]  x_hat  $\hat{x}$
 *
 * @return     $x$
 */
inline double Transform1D::forwardTransform(double x_hat)
{
	return x_nodes*gl.evalGLL(x_hat).leftCols<1>();
}

/**
 * @brief      Calculates the Jacobian of the transform at $\hat{x}$.
 *
 * @param[in]  x_hat  \hat{x}$
 *
 * @return     $J$
 */
inline double Transform1D::jacobian(double x_hat)
{
	return x_nodes*gl.evalGLL(x_hat).rightCols<1>();
}

#endif
