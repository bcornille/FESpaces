#include "Eigen/Core"
#include "Basis1D.hpp"

#ifndef _ElementTransforms_hpp
#define _ElementTransforms_hpp

using namespace Eigen;

// class Transform1D_Linear
// {
// 	public:
// 		Transform1D_Linear(double a = -1.0, double b = 1.0);
// 		~Transform1D_Linear() = default;
// 		double forwardTransform(double x_hat);
// 		double jacobian();
// 	private:
// 		const RowVector2d x_range;
// };

// Transform1D_Linear::Transform1D_Linear(double a, double b) : x_range({a, b}) {}

// inline double Transform1D_Linear::forwardTransform(double x_hat)
// {
// 	return ((x_range[1] - x_range[0])*x_hat + x_range[0] + x_range[1])/2.0;
// }

// inline double Transform1D_Linear::jacobian()
// {
// 	return (x_range[1] - x_range[0])/2.0;
// }

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

Transform1D::Transform1D(double a, double b) : x_nodes(2), gl(2)
{
	x_nodes[0] = a;
	x_nodes[1] = b;
}

Transform1D::Transform1D(RowVectorXd nodes) : x_nodes(nodes), gl(nodes.size()) {}

inline double Transform1D::forwardTransform(double x_hat)
{
	return x_nodes*gl.evalGLL(x_hat).leftCols<1>();
}

inline double Transform1D::jacobian(double x_hat)
{
	return x_nodes*gl.evalGLL(x_hat).rightCols<1>();
}

#endif
