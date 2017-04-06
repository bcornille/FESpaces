#include "Eigen/Core"

#ifndef _ElementTransforms_hpp
#define _ElementTransforms_hpp

using namespace Eigen;

class Trasform1D_Linear
{
	public:
		Trasform1D_Linear(double a = -1.0, double b = 1.0);
		~Trasform1D_Linear() = default;
		double forwardTrasform(double x_hat);
		double jacobian();
	private:
		RowVector2d x_range;
};

Trasform1D_Linear::Trasform1D_Linear(double a, dobule b) : x_range({a, b}) {}

inline Trasform1D_Linear::forwardTrasform(double x_hat)
{
	return ((b - a)*x_hat + a + b)/2.0;
}

inline Trasform1D_Linear::jacobian()
{
	return (b - a)/2.0;
}

#endif
