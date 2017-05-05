#include "json.hpp"
#include "Eigen/Core"
#include <vector>
#include <cmath>
#include "ElementTransforms.hpp"

#ifndef _Mesh1D_hpp
#define _Mesh1D_hpp

using namespace nlohmann;

class Mesh1D
{
	public:
		Mesh1D(json params);
		~Mesh1D() = default;
		// Transform1D_Linear getLinearTransform(int i);
		Transform1D getTransform(int i);
		std::vector<double> nodes();
	private:
		const int N_node;
		const int N_segments;
		std::vector<double> x;
		std::vector<RowVectorXd> inner_nodes;
};

Mesh1D::Mesh1D(nlohmann::json params) :
	N_node((int)params["N"]), N_segments((int)params["N"] + 1),
	x((int)params["N"] + 2), inner_nodes((int)params["N"] + 1)
{
	double x_len = (double)params["x_max"] - (double)params["x_min"];
	double deltaxi = 2.0/N_segments;
	double c = (double)params["c"];
	x[0] = (double)params["x_min"];
	for (int i = 1; i < N_segments; ++i)
	{
		double xi = -1.0 + i*deltaxi;
		x[i] = (1.0 + xi + c*sin(pi()*xi))*x_len/2.0;
	}
	x[N_segments] = (double)params["x_max"];
	GaussLobatto gl((int)params["map_order"] + 1);
	RowVectorXd temp_nodes(gl.getN());
	for (int i = 0; i < N_segments; ++i)
	{
		Transform1D segment(-1.0 + i*deltaxi, -1.0 + (i + 1)*deltaxi);
		for (int j = 0; j < gl.getN(); ++j)
		{
			double xi = segment.forwardTransform(gl.getNode(j));
			temp_nodes[j] = x_len*(1.0 + xi + c*sin(pi()*xi))/2.0;
		}
		inner_nodes[i] = temp_nodes;
	}
}

// Transform1D_Linear Mesh1D::getLinearTransform(int i)
// {
// 	return Transform1D_Linear(x[i], x[i+1]);
// }

Transform1D Mesh1D::getTransform(int i)
{
	return Transform1D(inner_nodes[i]);
}

inline std::vector<double> Mesh1D::nodes()
{
	return x;
}

#endif
