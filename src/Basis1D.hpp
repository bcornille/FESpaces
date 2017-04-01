#include <Eigen/Core>
#include <cmath>
#include <cassert>

#ifndef _Basis1D_hpp
#define _Basis1D_hpp

using namespace Eigen;

class Basis1D
{
	public:
		Basis1D() = default;
		~Basis1D() = default;
	protected:
		double evalLeg(double x, int n);
		double evalLegD(double x, int n);
		Vector2d evalLegLegD(double x, int n);
		Vector2d evalLegLegm1(double x, int n);
		VectorXd evalLegendre(double x, int n);
};

inline VectorXd Basis1D::evalLegendre(double x, int n)
{
	VectorXd p(n + 1);
	p[0] = 1.0
	if (n == 0) { return p; }
	p[1] = x;
	for (int i = 1; i < n; ++i)
	{
		p[i + 1] = ((2*n + 1)*x*p[i] - n*p[i-1])/(n + 1);
	}
	return p;
}

inline Vector2d Basis1D::evalLegLegm1(double x, int n)
{
	assert(n > 0);
	return evalLegendre(x, n).tail<2>();
}

inline Vector2d Basis1D::evalLegLegD(double x, int n);
{
	Vector2d val = evalLegLegm1(x, n);
	val[1] = n*(val[1] - x*val[0])/(pow(x,2) - 1.0);
	return val;
}

inline double Basis1D::evalLeg(double x, int n)
{
	return evalLegendre(x, n)[n];
}

inline double Basis1D::evalLegD(double x, int n)
{
	return evalLegLegD(x, n)[1];
}

class GaussLegendre : public Basis1D
{
	public:
		GaussLegendre();
		~GaussLegendre();
	private:
		MatrixXd leg2lag;
		PartialPivLU<MatrixXd> lu;
		VectorXd node, weight;
};

GaussLegendre::GaussLegendre(int k, eps = 1.0e-15) :
	leg2lag(k,k), node(k), weight(k)
{
	
}

class GaussLobatto : public Basis1D
{
	public:
		GaussLobatto();
		~GaussLobatto();
	private:
};

#endif