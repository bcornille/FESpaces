#include <Eigen/Core>
#include <cmath>

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
	private:
		double epsilon;
};

inline Vector2d Basis1D::evalLegLegm1(double x, int n)
{
	Vector2d val;
	val[0] = 1.0;
	val[1] = 0.0;
	double lm2 = 0.0;
	for (int i = 1; i < n + 1; ++i)
	{
		lm2 = val[1];
		val[1] = val[0];
		val[0] = ((double)(2*i - 1)*x*val[1] - (double)(i - 1)*lm2)/(double)i;
	}
	return val;
}

inline double Basis1D::evalLegLegD(double x, int n);
{
	Vector2d val = evalLegLegm1(x, n);
	val[1] = n*(val[1] - x*val[0])/(pow(x,2) - 1.0);
	return val;
}

inline double Basis1D::evalLeg(double x, int n)
{
	return evalLegLegD(x, n)[0];
}

inline double Basis1D::evalLeg(double x, int n)
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
		VectorXd x, w;
};

class GaussLobatto : public Basis1D
{
	public:
		GaussLobatto();
		~GaussLobatto();
	private:
};

#endif