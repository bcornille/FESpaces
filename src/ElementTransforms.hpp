#include "Eigen/Core"

#ifndef _ElementTransforms_hpp
#define _ElementTransforms_hpp

using namespace Eigen;

class Transform1D_Linear
{
	public:
		Transform1D_Linear(double a = -1.0, double b = 1.0);
		~Transform1D_Linear() = default;
		double forwardTransform(double x_hat);
		double jacobian();
	private:
		const RowVector2d x_range;
};

Transform1D_Linear::Transform1D_Linear(double a, double b) : x_range({a, b}) {}

inline double Transform1D_Linear::forwardTransform(double x_hat)
{
	return ((x_range[1] - x_range[0])*x_hat + x_range[0] + x_range[1])/2.0;
}

inline double Transform1D_Linear::jacobian()
{
	return (x_range[1] - x_range[0])/2.0;
}

#endif
