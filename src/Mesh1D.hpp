#include "json.hpp"
#include <vector>
#include "ElementTransforms.hpp"

#ifndef _Mesh1D_hpp
#define _Mesh1D_hpp

using namespace nlohmann;

class Mesh1D
{
	public:
		Mesh1D(json params);
		~Mesh1D() = default;
		Transform1D_Linear getLinearTransform(int i);
		std::vector<double> nodes();
	private:
		std::vector<double> x;
		const int N_segments;
		const int N_node;
};

Mesh1D::Mesh1D(nlohmann::json params) :
	N_node((int)params["N"]), N_segments((int)params["N"] + 1), x((int)params["N"] + 2)
{
	double deltax = (((double)params["x_max"] - (double)params["x_min"])
		/(x.size() - 1));
	for (int i = 0; i < x.size(); ++i)
	{
		x[i] = (int)params["x_min"] + i*deltax;
	}
}

Transform1D_Linear Mesh1D::getLinearTransform(int i)
{
	return Transform1D_Linear(x[i], x[i+1]);
}

inline std::vector<double> Mesh1D::nodes()
{
	return x;
}

#endif
