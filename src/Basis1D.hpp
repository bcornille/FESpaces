#include "Eigen/Core"
#include "Eigen/Dense"
#include <cmath>

#ifndef _Basis1D_hpp
#define _Basis1D_hpp

using namespace Eigen;

constexpr double pi() { return std::atan(1)*4; }

class Basis1D
{
	public:
		Basis1D() = default;
		virtual ~Basis1D() = default;
	protected:
		// double evalLeg(double x, int n);
		// double evalLegD(double x, int n);
		RowVector2d evalLegLegD(double x, int n);
		Vector2d evalLegm1Leg(double x, int n);
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
		p(i+1, 1) = (2*i + 1)*p(i, 0) + p(i-1, 1);
	}
	return p;
}

inline VectorXd Basis1D::evalLegendre(double x, int n)
{
	return evalLegendreD(x, n).leftCols<1>();
}

inline Vector2d Basis1D::evalLegm1Leg(double x, int n)
{
	return evalLegendreD(x, n).bottomLeftCorner<2,1>();
}

inline RowVector2d Basis1D::evalLegLegD(double x, int n)
{
	return evalLegendreD(x, n).bottomRows<1>();
}

// inline double Basis1D::evalLeg(double x, int n)
// {
// 	return evalLegendreD(x, n)(n,0);
// }

// inline double Basis1D::evalLegD(double x, int n)
// {
// 	return evalLegendreD(x, n)(n,1);
// }

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
		double x = cos((2.0*i - 1.0)*pi()/(2.0*k));
		do
		{
			double x_old = x;
			p = evalLegLegD(x, k);
			x = x_old - p[0]/p[1];
			err = std::abs(x_old - x);
		} while (err > eps);
		node[i-1] = -x;
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
	return lu.solve(evalLegendreD(x, n_nodes - 1));
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

class GaussLobatto : public Basis1D
{
	public:
		GaussLobatto(int k = 1, double eps = 1.0e-15);
		virtual ~GaussLobatto() = default;
		MatrixX2d evalGLL(double x);
		double getNode(int i);
		double getWeight(int i);
		int getN();
	private:
		const int n_nodes;
		MatrixXd leg2lag;
		PartialPivLU<MatrixXd> lu;
		VectorXd node, weight;
};

GaussLobatto::GaussLobatto(int k, double eps) :
	n_nodes(k), leg2lag(k,k), lu(k), node(k), weight(k)
{
	GaussLegendre gl(k-1);
	node[0] = -1.0;
	node[k-1] = 1.0;
	weight[0] = 2.0/(k*(k-1));
	weight[k-1] = weight[0];
	for (int i = 2; i <= (k + 1)/2; ++i)
	{
		double err;
		Vector2d p;
		double x = 0.5*(gl.getNode(i-1) + gl.getNode(i-2));
		do
		{
			double x_old = x;
			p = evalLegm1Leg(x, k - 1);
			double g = k*(x*p[1] - p[0]);
			double dg = k*(k + 1)*p[1];
			x = x_old - g/dg;
			err = std::abs(x_old - x);
		} while (err > eps);
		node[i-1] = x;
		node[k-i] = -x;
		weight[i-1] = weight[0]/pow(p[1],2);
		weight[k-i] = weight[i-1];
	}
	for (int i = 0; i < k; ++i)
	{
		leg2lag.col(i) = evalLegendre(node[i], k - 1);
	}
	lu.compute(leg2lag);
}

inline MatrixX2d GaussLobatto::evalGLL(double x)
{
	return lu.solve(evalLegendreD(x, n_nodes - 1));
}

inline double GaussLobatto::getNode(int i)
{
	return node[i];
}

inline double GaussLobatto::getWeight(int i)
{
	return weight[i];
}

inline int GaussLobatto::getN()
{
	return n_nodes;
}

class EdgeFunction : private GaussLobatto
{
	public:
		EdgeFunction(int k = 1, double eps = 1e-15);
		~EdgeFunction() = default;
		VectorXd evalEF(double x);
		int getN();
	private:
		const int n_segments;
};

EdgeFunction::EdgeFunction(int k, double eps) :
	GaussLobatto(k + 1, eps), n_segments(k) {}

inline VectorXd EdgeFunction::evalEF(double x)
{
	VectorXd ef(n_segments);
	MatrixX2d gl = evalGLL(x);
	ef[0] = -gl(0, 1);
	for (int i = 1; i < n_segments; ++i)
	{
		ef[i] = ef[i-1] - gl(i, 1);
	}
	return ef;
}

inline int EdgeFunction::getN()
{
	return n_segments;
}

#endif
