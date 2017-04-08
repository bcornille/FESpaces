#include <cmath>

#ifndef _ForcingFunctions_hpp
#define _ForcingFunctions_hpp

class Force1D
{
	public:
		// Force1D() = default;
		virtual ~Force1D() = default;
		virtual double f(double x) = 0;
};

class ExpX3pX : public Force1D
{
	public:
		// ExpX1mX() = default;
		~ExpX3pX() = default;
		double f(double x);
};

double ExpX3pX::f(double x)
{
	return exp(x)*x*(x + 3.0);
}

class Two : public Force1D
{
	public:
		// Two() = default;
		~Two() = default;
		double f(double x);
};

double Two::f(double x)
{
	return 2.0;
}

#endif
