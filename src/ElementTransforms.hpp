#include "Eigen/Core"
#include "Eigen/LU"
#include "ReferenceElements.hpp"

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
		GaussLobatto gll;
};

/**
 * @brief      Constructs the object. Creates a linear transformation.
 *
 * @param[in]  a     Minimum $x$ value of the segment.
 * @param[in]  b     Maximum $x$ value of the segment.
 */
Transform1D::Transform1D(double a, double b) : x_nodes(2), gll(2)
{
	x_nodes[0] = a;
	x_nodes[1] = b;
}

/**
 * @brief      Constructs the object. Creates an arbitrary transformation.
 *
 * @param[in]  nodes  The nodes corresponding to the transformed GLL points.
 */
Transform1D::Transform1D(RowVectorXd nodes) : x_nodes(nodes), gll(nodes.size()) {}

/**
 * @brief      Calculates the forward transformation $x=F(\hat{x})$.
 *
 * @param[in]  x_hat  $\hat{x}$
 *
 * @return     $x$
 */
inline double Transform1D::forwardTransform(double x_hat)
{
	return x_nodes*gll.evalGLL(x_hat).leftCols<1>();
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
	return x_nodes*gll.evalGLL(x_hat).rightCols<1>();
}

class Transform2D
{
	public:
		Transform2D(Vector2d x_min = Vector2d::Constant(-1.0),
			Vector2d x_max = Vector2d::Constant(1.0));
		Transform2D(Matrix2Xd nodes, int order);
		~Transform2D() = default;
		Vector2d forwardTransform(Vector2d x_hat);
		Matrix2d jacobianMatrix(Vector2d x_hat);
		Matrix2d jacobianInverse(Vector2d x_hat);
		double jacobian(Vector2d x_hat);
	private:
		int n_map;
		Matrix2Xd x_nodes;
		H1_2D shape;
};

Transform2D::Transform2D(Vector2d x_min, Vector2d x_max) :
	n_map(4), x_nodes(2, 4), shape(1)
{
	x_nodes.col(0) = x_min;
	x_nodes.col(1) << x_max[0], x_min[1];
	x_nodes.col(2) << x_min[0], x_max[1];
	x_nodes.col(3) = x_max;
}

Transform2D::Transform2D(Matrix2Xd nodes, int order) :
	n_map((order + 1)*(order + 1)), x_nodes(nodes), shape(order) {}

inline Vector2d Transform2D::forwardTransform(Vector2d x_hat)
{
	return x_nodes*shape.eval(x_hat);
}

inline Matrix2d Transform2D::jacobianMatrix(Vector2d x_hat)
{
	return x_nodes*shape.evalGrad(x_hat);
}

inline Matrix2d Transform2D::jacobianInverse(Vector2d x_hat)
{
	return (x_nodes*shape.evalGrad(x_hat)).inverse();
}

inline double Transform2D::jacobian(Vector2d x_hat)
{
	return jacobianMatrix(x_hat).determinant();
}

#endif
