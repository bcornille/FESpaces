#define _USE_MATH_DEFINES
#include "Eigen/Core"
#include "Eigen/Dense"
#include <cmath>

#ifndef _Basis1D_hpp
#define _Basis1D_hpp

using namespace Eigen;

class Basis1D
{
	public:
		Basis1D() = default;
		virtual ~Basis1D() = default;
	protected:
		double evalLeg(double x, int n);
		double evalLegD(double x, int n);
		RowVector2d evalLegLegD(double x, int n);
		Vector2d evalLegLegm1(double x, int n);
		VectorXd evalLegendre(double x, int n);
		MatrixX2d evalLegendreD(double x, int n);
};

inline MatrixX2d Basis1D::evalLegendreD(double x, int n)
{
	MatrixX2d p(n + 1, 2);
	p(0, 0) = 1.0;
	p(0, 1) = 0.0;
	if (n == 0) { return p; }
	p(1, 0) = x;
	p(1, 1) = 1.0;
	for (int i = 1; i < n; ++i)
	{
		p(i+1, 0) = ((2*i + 1)*x*p(i, 0) - i*p(i-1, 0))/(i + 1);
		p(i+1, 1) = (4*i + 2)*p(i, 0) + p(i-1, 1);
	}
	return p;
}

inline VectorXd Basis1D::evalLegendre(double x, int n)
{
	return evalLegendreD(x, n).leftCols<1>();
}

inline Vector2d Basis1D::evalLegLegm1(double x, int n)
{
	return evalLegendreD(x, n).bottomLeftCorner<2,1>();
}

inline RowVector2d Basis1D::evalLegLegD(double x, int n)
{
	return evalLegendreD(x, n).bottomRows<1>();
}

inline double Basis1D::evalLeg(double x, int n)
{
	return evalLegendreD(x, n)(n,0);
}

inline double Basis1D::evalLegD(double x, int n)
{
	return evalLegendreD(x, n)(n,1);
}

class GaussLegendre : public Basis1D
{
	public:
		GaussLegendre(int k = 1, double eps = 1.0e-15);
		~GaussLegendre() = default;
		MatrixX2d evalGL(double x);
		double getNode(int i);
		double getWeight(int i);
		int getN();
	private:
		const int n_nodes;
		MatrixXd leg2lag;
		PartialPivLU<MatrixXd> lu;
		VectorXd node, weight;
};

GaussLegendre::GaussLegendre(int k, double eps) :
	n_nodes(k), leg2lag(k,k), lu(k), node(k), weight(k)
{
	for (int i = 1; i <= (k + 1)/2; ++i)
	{
		double err;
		RowVector2d p;
		double x = cos((2.0*i - 1.0)*M_PI/(2.0*k));
		do
		{
			double x_old = x;
			p = evalLegLegD(x, k);
			x = x_old - p[0]/p[1];
			err = std::abs(x_old - x);
		} while (err > eps);
		node[i-1] = x;
		node[k-i] = x;
		weight[i-1] = 2.0/((1.0 - pow(x,2))*pow(p[1],2));
		weight[k-i] = weight[i-1];
	}
	for (int i = 0; i < k; ++i)
	{
		leg2lag.col(i) = evalLegendre(node[i], k - 1);
	}
	lu.compute(leg2lag);
}

inline MatrixX2d GaussLegendre::evalGL(double x)
{
	return lu.solve(evalLegendreD(x, n_nodes));
}

inline double GaussLegendre::getNode(int i)
{
	return node[i];
}

inline double GaussLegendre::getWeight(int i)
{
	return weight[i];
}

inline int GaussLegendre::getN()
{
	return n_nodes;
}

// class GaussLobatto : public Basis1D
// {
// 	public:
// 		GaussLobatto() = default;
// 		~GaussLobatto() = default;
// 	private:
// };

#endif