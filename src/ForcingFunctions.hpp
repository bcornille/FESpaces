#include <cmath>

#ifndef _ForcingFunctions_hpp
#define _ForcingFunctions_hpp

class Force1D
{
	public:
		Force1D() = default;
		virtual ~Force1D() = default;
		virtual double f(double x);
};

class ExpX1mX : Force1D
{
	public:
		ExpX1mX() = default;
		~ExpX1mX() = default;
		double f(double x);
};

double ExpX1mX::f(double x)
{
	return exp(x)*x*(1.0 - x);
}

class Two : Force1D
{
	public:
		Two() = default;
		~Two() = default;
		double f(double x);
};

double Two::f(double x)
{
	return 2.0;
}

#endif
