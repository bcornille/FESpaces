#include "Eigen/Core"
#include "Eigen/Dense"
#include <cmath>

#ifndef _Basis1D_hpp
#define _Basis1D_hpp

using namespace Eigen;

/**
 * @brief      Gives the value of pi. constexpr means it is evaluated at compile
 *             time.
 *
 * @return     pi = std::atan(1)*4
 */
constexpr double pi() { return std::atan(1)*4; }

/**
 * @brief      Class for developing polynomial bases in 1D. Provides functions
 *             for evaluating the Legendre polynomials, which can be used to
 *             construct other bases.
 */
class Basis1D
{
	public:
		Basis1D() = default;
		virtual ~Basis1D() = default;
	protected:
		RowVector2d evalLegLegD(double x, int n);
		Vector2d evalLegm1Leg(double x, int n);
		VectorXd evalLegendre(double x, int n);
		MatrixX2d evalLegendreD(double x, int n);
};

/**
 * @brief      Evaluate the Legendre polynomials and their derivatives up to a
 *             certain degree.
 *
 * @param[in]  x     Coordinate value `-1 < x < 1`.
 * @param[in]  n     Order of highest Legendre polynomial.
 *
 * @return    `n+1`-by-`2` matrix with the first column holding the values of the
 *             Legendre polynomials and the second column holding the valuse of
 *             the derivatives of the Legendre polynomials.
 */
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

/**
 * @brief      Evaluate the Legendre polynomials up to a certain degree.
 *
 * @param[in]  x     Coordinate value `-1 < x < 1`.
 * @param[in]  n     Order of highest Legendre polynomial.
 *
 * @return     Vector of length `n+1` holding the values of the Legendre
 *             polynomials.
 */
inline VectorXd Basis1D::evalLegendre(double x, int n)
{
	return evalLegendreD(x, n).leftCols<1>();
}

/**
 * @brief      Evaluate a Legendre polynomial of a certain degre and degree
 *             minus 1.
 *
 * @param[in]  x     Coordinate value `-1 < x < 1`.
 * @param[in]  n     Order of highest Legendre polynomial.
 *
 * @return     Vector of length `2` containing the two values.
 */
inline Vector2d Basis1D::evalLegm1Leg(double x, int n)
{
	return evalLegendreD(x, n).bottomLeftCorner<2,1>();
}

/**
 * @brief      Evaluate a Legendre polynomial of a certain degree and its
 *             derivative.
 *
 * @param[in]  x     Coordinate value `-1 < x < 1`.
 * @param[in]  n     Order of Legendre polynomial.
 *
 * @return     Row vector of length `2` containing the two values.
 */
inline RowVector2d Basis1D::evalLegLegD(double x, int n)
{
	return evalLegendreD(x, n).bottomRows<1>();
}

/**
 * @brief      Class for the Lagrange cardinal basis at the Gauss-Legendre (GL)
 *             points.
 */
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
		int n_nodes;
		MatrixXd leg2lag;
		PartialPivLU<MatrixXd> lu;
		VectorXd node, weight;
};

/**
 * @brief      Constructs the object. During construction the GL nodes and
 *             weights are calculated. The values of the Legendre polynomials up
 *             to degree `k - 1` at each node are also stored in a matrix. To
 *             evaluate the cardinal functions, this matrix is inverted on the
 *             vector containing the evalutation of the Legendre polynomials at
 *             the desired point.
 *
 * @param[in]  k     The number of GL points to use.
 * @param[in]  eps   The required accuracy of the GL points.
 */
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

/**
 * @brief      Evaluates the Lagrange GL cardinal fuctions and their derivatives
 *             at a desired location.
 *
 * @param[in]  x     Coordinate value `-1 < x < 1`.
 *
 * @return     `k`-by-`2` matrix of the cardinal function values and their
 *             derivatives.
 */
inline MatrixX2d GaussLegendre::evalGL(double x)
{
	return lu.solve(evalLegendreD(x, n_nodes - 1));
}

/**
 * @brief      Gets the node.
 *
 * @param[in]  i     Index of desired node.
 *
 * @return     The node location.
 */
inline double GaussLegendre::getNode(int i)
{
	return node[i];
}

/**
 * @brief      Gets the weight.
 *
 * @param[in]  i     Index of desired node.
 *
 * @return     The weight.
 */
inline double GaussLegendre::getWeight(int i)
{
	return weight[i];
}

/**
 * @brief      Gets the number of nodes.
 *
 * @return     The number of nodes.
 */
inline int GaussLegendre::getN()
{
	return n_nodes;
}

/**
 * @brief      Class for the Lagrange cardinal basis at the
 *             Gauss-Legendre-Lobatto (GLL) points.
 */
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
		int n_nodes;
		MatrixXd leg2lag;
		PartialPivLU<MatrixXd> lu;
		VectorXd node, weight;
};

/**
 * @brief      Constructs the object. During construction the GLL nodes and
 *             weights are calculated. The values of the Legendre polynomials up
 *             to degree `k - 1` at each node are also stored in a matrix. To
 *             evaluate the cardinal functions, this matrix is inverted on the
 *             vector containing the evalutation of the Legendre polynomials at
 *             the desired point.
 *
 * @param[in]  k     The number of GLL points to use.
 * @param[in]  eps   The required accuracy of the GLL points.
 */
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

/**
 * @brief      Evaluates the Lagrange GLL cardinal fuctions and their
 *             derivatives at a desired location.
 *
 * @param[in]  x     Coordinate value `-1 < x < 1`.
 *
 * @return     `k`-by-`2` matrix of the cardinal function values and their
 *             derivatives.
 */
inline MatrixX2d GaussLobatto::evalGLL(double x)
{
	return lu.solve(evalLegendreD(x, n_nodes - 1));
}

/**
 * @brief      Gets the node.
 *
 * @param[in]  i     Index of desired node.
 *
 * @return     The node.
 */
inline double GaussLobatto::getNode(int i)
{
	return node[i];
}

/**
 * @brief      Gets the weight.
 *
 * @param[in]  i     Index of desire node
 *
 * @return     The weight.
 */
inline double GaussLobatto::getWeight(int i)
{
	return weight[i];
}

/**
 * @brief      Gets the number of nodes.
 *
 * @return     The number of nodes.
 */
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
		int n_segments;
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
