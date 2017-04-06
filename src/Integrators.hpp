#include "json.hpp"
#include "ReferenceElements.hpp"
#include "ElementTransforms.hpp"

#ifndef _Integrators_hpp
#define _Integrators_hpp

using namespace nlohmann;
using namespace Eigen;

class Integrator1D
{
	public:
		Integrator1D(json params);
		~Integrator1D() = default;
	private:
		GaussLegendre gl;
};

Integrator::Integrator(nlohmann::json params) : gl((int)params["N"]) {}

#endif