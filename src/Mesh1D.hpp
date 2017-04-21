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
	double deltax = (((double)params["x_max"] - (double)params["x_min"])
		/(x.size() - 1));
	for (int i = 0; i < x.size(); ++i)
	{
		x[i] = (int)params["x_min"] + i*deltax;
	}
	GaussLobatto gl((int)params["map_order"] + 1);
	RowVectorXd temp_nodes(gl.getN());
	double c = (double)params["c"];
	for (int i = 0; i < inner_nodes.size(); ++i)
	{
		for (int j = 0; j < temp_nodes.size(); ++j)
		{
			double xi = gl.getNode(j);
			temp_nodes[j] = deltax*(1.0 + xi + c*cos(pi()/2.0*xi+pi()*(i+1)))/2.0;
		}
		inner_nodes[i] = temp_nodes + RowVectorXd::Constant(temp_nodes.size(), x[i]);
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
