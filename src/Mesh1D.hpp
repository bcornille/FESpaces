#include "json.hpp"
#include "Eigen/Core"
#include <vector>
#include <cmath>
#include "ElementTransforms.hpp"

#ifndef _Mesh1D_hpp
#define _Mesh1D_hpp

using namespace nlohmann;

/**
 * @brief      Class for one dimensional mesh. Can include some curvature.
 */
class Mesh1D
{
	public:
		Mesh1D(json params);
		~Mesh1D() = default;
		Transform1D getTransform(int i);
		std::vector<double> nodes();
	private:
		const int N_node;
		const int N_segments;
		std::vector<double> x;
		std::vector<RowVectorXd> inner_nodes;
};

/**
 * @brief      Constructs the object. Produces `"N" + 1` elements. This mesh is
 *             distorted by the mapping from $\xi \in [-1, 1]$ to $1 + \xi +
 *             c\sin(\pi \xi)$.
 *
 * @param[in]  params  The parameters of "Mesh". Required to have fields
 *                     "x_min", "x_max", "N", "c", and "map_order".
 */
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

/**
 * @brief      Gets the one dimensional transform for a desired element.
 *
 * @param[in]  i     Element index.
 *
 * @return     The transform.
 */
Transform1D Mesh1D::getTransform(int i)
{
	return Transform1D(inner_nodes[i]);
}

/**
 * @brief      Gets the set of "nodes" for the mesh.
 *
 * @return     The nodes.
 */
inline std::vector<double> Mesh1D::nodes()
{
	return x;
}

#endif
