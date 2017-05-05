#include "Basis1D.hpp"

#ifndef _ForcingFunctions_hpp
#define _ForcingFunctions_hpp

using namespace Eigen;

class Force1D
{
	public:
		// Force1D() = default;
		virtual ~Force1D() = default;
		virtual double f(double x) = 0;
		virtual double sol(double x) = 0;
};

class ExpX3pX : public Force1D
{
	public:
		// ExpX1mX() = default;
		~ExpX3pX() = default;
		double f(double x);
		double sol(double x);
};

double ExpX3pX::f(double x)
{
	return exp(x)*x*(x + 3.0);
}

double ExpX3pX::sol(double x)
{
	return exp(x)*x*(1.0 - x);
}

class Two : public Force1D
{
	public:
		// Two() = default;
		~Two() = default;
		double f(double x);
		double sol(double x);
};

double Two::f(double x)
{
	return 2.0;
}

double Two::sol(double x)
{
	return x*(1.0 - x);
}

class Force2D
{
	public:
		// Force2D() = default;
		virtual ~Force2D() = default;
		virtual double f(Vector2d x) = 0;
		virtual double sol(Vector2d x) = 0;
};

class SinXSinY : public Force2D
{
	public:
		~SinXSinY() = default;
		double f(Vector2d x);
		double sol(Vector2d x);
};

double SinXSinY::f(Vector2d x)
{
	return 2*pow(pi(), 2)*sin(pi()*x[0])*sin(pi()*x[1]);
}

double SinXSinY::sol(Vector2d x)
{
	return sin(pi()*x[0])*sin(pi()*x[1]);
}

#endif
