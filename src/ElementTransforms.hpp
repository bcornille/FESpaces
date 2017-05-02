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
		GaussLobatto gll;
};

Transform1D::Transform1D(double a, double b) : x_nodes(2), gll(2)
{
	x_nodes[0] = a;
	x_nodes[1] = b;
}

Transform1D::Transform1D(RowVectorXd nodes) : x_nodes(nodes), gll(nodes.size()) {}

inline double Transform1D::forwardTransform(double x_hat)
{
	return x_nodes*gll.evalGLL(x_hat).leftCols<1>();
}

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
	private:
		int n_map;
		Matrix2Xd x_nodes;
		GaussLobatto gll;
};

Transform2D::Transform2D(Vector2d x_min, Vector2d x_max) :
	n_map(4), x_nodes(4, 2), gll(2)
{
	x_nodes.col(0) = x_min;
	x_nodes.col(1) << x_max[0], x_min[1];
	x_nodes.col(2) << x_min[0], x_max[1];
	x_nodes.col(3) = x_max;
}

Transform2D::Transform2D(Matrix2Xd nodes, int order) :
	n_map((order + 1)*(order + 1)), x_nodes(nodes), gll(order + 1) {}

inline Vector2d Transform2D::forwardTransform(Vector2d x_hat)
{
	Vector2d trans;
	MatrixXd mat_vals = (gll.evalGLL(x_hat[0]).leftCols<1>()
		*gll.evalGLL(x_hat[1]).leftCols<1>().transpose());
	VectorXd vec_vals = Map<VectorXd>(mat_vals.data(), n_map);
	trans = x_nodes*vec_vals;
	return trans;
}

#endif
