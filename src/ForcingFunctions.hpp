#include <cmath>

#ifndef _ForcingFunctions_hpp
#define _ForcingFunctions_hpp

/**
 * @brief      Class for one dimensional forcing fuction. This is a pure virtual
 *             class to serve as an interface for possible forcing functions.
 */
class Force1D
{
	public:
		// Force1D() = default;
		virtual ~Force1D() = default;
		virtual double f(double x) = 0;
		virtual double sol(double x) = 0;
};

/**
 * @brief      Class for the forcing function $f(x) = e^x x (3 + x)$, which has
 *             the solution $p(x)=e^x x (1 - x)$.
 */
class ExpX3pX : public Force1D
{
	public:
		// ExpX1mX() = default;
		~ExpX3pX() = default;
		double f(double x);
		double sol(double x);
};

/**
 * @brief      $f(x) = e^x x (3 + x)$
 *
 * @param[in]  x     $x$
 *
 * @return     f
 */
double ExpX3pX::f(double x)
{
	return exp(x)*x*(x + 3.0);
}

/**
 * @brief      $p(x) = e^x x (1 - x)$
 *
 * @param[in]  x     $x$
 *
 * @return     $p$
 */
double ExpX3pX::sol(double x)
{
	return exp(x)*x*(1.0 - x);
}

/**
 * @brief      Class for the forcing function $f(x) = 2$, which has the solution
 *             $p(x) = x (1 - x)$.
 */
class Two : public Force1D
{
	public:
		// Two() = default;
		~Two() = default;
		double f(double x);
		double sol(double x);
};

/**
 * @brief      $f(x) = 2$
 *
 * @param[in]  x     $x$
 *
 * @return     $f$
 */
double Two::f(double x)
{
	return 2.0;
}

/**
 * @brief      $p(x) = x (1 - x)$
 *
 * @param[in]  x     $x$
 *
 * @return     $p$
 */
double Two::sol(double x)
{
	return x*(1.0 - x);
}

#endif
