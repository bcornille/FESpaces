#include "json.hpp"
#include "ReferenceElements.hpp"
#include "ElementTransforms.hpp"
#include "ForcingFunctions.hpp"

#ifndef _Integrators_hpp
#define _Integrators_hpp

using namespace nlohmann;
using namespace Eigen;

/**
 * @brief      Class for one dimensional integrators. Gauss-Legendre quadradure
 *             is used for integration.
 */
class Integrator1D
{
	public:
		Integrator1D(json params);
		~Integrator1D() = default;
		MatrixXd mass(H1_1D u, H1_1D v, Transform1D t);
		MatrixXd mass(L2_1D u, L2_1D v, Transform1D t);
		MatrixXd mass(L2_1D_EF u, L2_1D_EF v, Transform1D t);
		MatrixXd grad(H1_1D u, L2_1D v, Transform1D t);
		MatrixXd grad(L2_1D u, H1_1D v, Transform1D t);
		MatrixXd grad(H1_1D u, L2_1D_EF v, Transform1D t);
		MatrixXd grad(L2_1D_EF u, H1_1D v, Transform1D t);
		MatrixXd laplace(H1_1D u, H1_1D v, Transform1D t);
		VectorXd force(const std::shared_ptr<Force1D>& f, H1_1D v, Transform1D t);
		VectorXd force(const std::shared_ptr<Force1D>& f, L2_1D v, Transform1D t);
		VectorXd force(const std::shared_ptr<Force1D>& f, L2_1D_EF v, Transform1D t);
		double error(const std::shared_ptr<Force1D>& f, H1_1D v,
			VectorXd coeffs, Transform1D t);
		double error(const std::shared_ptr<Force1D>& f, L2_1D v,
			VectorXd coeffs, Transform1D t);
		double error(const std::shared_ptr<Force1D>& f, L2_1D_EF v,
			VectorXd coeffs, Transform1D t);
	private:
		GaussLegendre gl;
};

/**
 * @brief      Constructs the object.
 *
 * @param[in]  params  The parameters of the "Integrator". Must have attribute "N".
 */
Integrator1D::Integrator1D(json params) : gl((int)params["N"]) {}

/**
 * @brief      Local mass matrix of $\left<u, v\right>$ where $u, v \in H^1$.
 *
 * @param[in]  u     $H^1$ element on $[-1, 1]$.
 * @param[in]  v     $H^1$ element on $[-1, 1]$.
 * @param[in]  t     1D transformation.
 *
 * @return     Local mass matrix.
 */
inline MatrixXd Integrator1D::mass(H1_1D u, H1_1D v, Transform1D t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd u_vals = u.eval(x_hat);
		VectorXd v_vals = v.eval(x_hat);
		matrix += v_vals*u_vals.transpose()*t.jacobian(x_hat)*gl.getWeight(i);
	}
	return matrix;
}

/**
 * @brief      Local mass matrix of $\left<u, v\right>$ where $u, v \in L^2$.
 *
 * @param[in]  u     $L^2$ element on $[-1, 1]$.
 * @param[in]  v     $L^2$ element on $[-1, 1]$.
 * @param[in]  t     1D transformation.
 *
 * @return     Local mass matrix.
 */
inline MatrixXd Integrator1D::mass(L2_1D u, L2_1D v, Transform1D t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd u_vals = u.eval(x_hat);
		VectorXd v_vals = v.eval(x_hat);
		matrix += v_vals*u_vals.transpose()*t.jacobian(x_hat)*gl.getWeight(i);
	}
	return matrix;
}

inline MatrixXd Integrator1D::mass(L2_1D_EF u, L2_1D_EF v, Transform1D t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd u_vals = u.eval(x_hat)/t.jacobian(x_hat);
		VectorXd v_vals = v.eval(x_hat)/t.jacobian(x_hat);
		matrix += v_vals*u_vals.transpose()*t.jacobian(x_hat)*gl.getWeight(i);
	}
	return matrix;
}

/**
 * @brief      Local integration of the bilinear operator $\left<u', v\right>$
 *             where $u \in H^1$ and $v \in L^2$.
 *
 * @param[in]  u     $H^1$ element on $[-1, 1]$
 * @param[in]  v     $L^2$ element on $[-1, 1]$
 * @param[in]  t     1D transformation.
 *
 * @return     Local matrix of the derivative operator.
 */
inline MatrixXd Integrator1D::grad(H1_1D u, L2_1D v, Transform1D t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd u_vals = u.evalD(x_hat);
		VectorXd v_vals = v.eval(x_hat);
		matrix += v_vals*u_vals.transpose()*gl.getWeight(i);
	}
	return matrix;
}

/**
 * @brief      Local integration of the bilinear operator $\left<u, v'\right>$
 *             where $u \in L^2$ and $v \in H^1$.
 *
 * @param[in]  u     $L^2$ element on $[-1, 1]$.
 * @param[in]  v     $H^1$ element on $[-1, 1]$.
 * @param[in]  t     1D transformation.
 *
 * @return     Local matrix of the derivative operator.
 */
inline MatrixXd Integrator1D::grad(L2_1D u, H1_1D v, Transform1D t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd u_vals = u.eval(x_hat);
		VectorXd v_vals = v.evalD(x_hat);
		matrix += v_vals*u_vals.transpose()*gl.getWeight(i);
	}
	return matrix;
}

inline MatrixXd Integrator1D::grad(H1_1D u, L2_1D_EF v, Transform1D t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd u_vals = u.evalD(x_hat);
		VectorXd v_vals = v.eval(x_hat)/t.jacobian(x_hat);
		matrix += v_vals*u_vals.transpose()*gl.getWeight(i);
	}
	return matrix;
}

inline MatrixXd Integrator1D::grad(L2_1D_EF u, H1_1D v, Transform1D t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd u_vals = u.eval(x_hat)/t.jacobian(x_hat);
		VectorXd v_vals = v.evalD(x_hat);
		matrix += v_vals*u_vals.transpose()*gl.getWeight(i);
	}
	return matrix;
}

/**
 * @brief      Local integration of the bilinear operator $\left<u', v'\right>$
 *             where $u, v \in H^1$.
 *
 * @param[in]  u     $H^1$ element on $[-1, 1]$.
 * @param[in]  v     $H^1$ element on $[-1, 1]$.
 * @param[in]  t     1D transformation.
 *
 * @return     Local mattrix of the double derivative operator.
 */
inline MatrixXd Integrator1D::laplace(H1_1D u, H1_1D v, Transform1D t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd u_vals = u.evalD(x_hat);
		VectorXd v_vals = v.evalD(x_hat);
		matrix += v_vals*u_vals.transpose()/t.jacobian(x_hat)*gl.getWeight(i);
	}
	return matrix;
}

/**
 * @brief      Local integration of the linear operator $\int fv$ where $v \in
 *             H^1$.
 *
 * @param[in]  f     Forcing function.
 * @param[in]  v     $H^1$ element on $[-1, 1]$.
 * @param[in]  t     1D transformation.
 *
 * @return     Local right hand side.
 */
inline VectorXd Integrator1D::force(const std::shared_ptr<Force1D>& f, H1_1D v,
	Transform1D t)
{
	VectorXd rhs(v.dofs());
	rhs.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		double x = t.forwardTransform(x_hat);
		VectorXd v_vals = v.eval(x_hat);
		rhs += f->f(x)*v_vals*t.jacobian(x_hat)*gl.getWeight(i);
	}
	return rhs;
}

/**
 * @brief      Local integration of the linear operator $\int fv$ where $v \in
 *             L^2$.
 *
 * @param[in]  f     Forcing function.
 * @param[in]  v     $L^2$ element on $[-1, 1]$.
 * @param[in]  t     1D transformation.
 *
 * @return     Local right hand side.
 */
inline VectorXd Integrator1D::force(const std::shared_ptr<Force1D>& f, L2_1D v,
	Transform1D t)
{
	VectorXd rhs(v.dofs());
	rhs.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		double x = t.forwardTransform(x_hat);
		VectorXd v_vals = v.eval(x_hat);
		rhs += f->f(x)*v_vals*t.jacobian(x_hat)*gl.getWeight(i);
	}
	return rhs;
}

inline VectorXd Integrator1D::force(const std::shared_ptr<Force1D>& f, L2_1D_EF v,
	Transform1D t)
{
	VectorXd rhs(v.dofs());
	rhs.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		double x = t.forwardTransform(x_hat);
		VectorXd v_vals = v.eval(x_hat)/t.jacobian(x_hat);
		rhs += f->f(x)*v_vals*t.jacobian(x_hat)*gl.getWeight(i);
	}
	return rhs;
}

/**
 * @brief      Local integration of the error $\left\|p - p_h\right\|^2_2$.
 *
 * @param[in]  f       Forcing function with solution.
 * @param[in]  v       $H^1$ element on $[-1, 1]$.
 * @param[in]  coeffs  The coefficients of the solution vector for the given
 *                     element.
 * @param[in]  t       1D transformation.
 *
 * @return     Local error squared.
 */
inline double Integrator1D::error(const std::shared_ptr<Force1D>& f, H1_1D v,
	VectorXd coeffs, Transform1D t)
{
	double err = 0.0;
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		double x = t.forwardTransform(x_hat);
		VectorXd v_vals = v.eval(x_hat);
		err += pow(f->sol(x) - v_vals.dot(coeffs), 2)*t.jacobian(x_hat)*gl.getWeight(i);
	}
	return err;
}

/**
 * @brief      Local integration of the error $\left\|p - p_h\right\|^2_2$.
 *
 * @param[in]  f       Forcing function with solution.
 * @param[in]  v       $L^2$ element on $[-1, 1]$.
 * @param[in]  coeffs  The coefficients of the solution vector for the given
 *                     element.
 * @param[in]  t       1D transformation.
 *
 * @return     Local error squared.
 */
inline double Integrator1D::error(const std::shared_ptr<Force1D>& f, L2_1D v,
	VectorXd coeffs, Transform1D t)
{
	double err = 0.0;
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		double x = t.forwardTransform(x_hat);
		VectorXd v_vals = v.eval(x_hat);
		err += pow(f->sol(x) - v_vals.dot(coeffs), 2)*t.jacobian(x_hat)*gl.getWeight(i);
	}
	return err;
}

inline double Integrator1D::error(const std::shared_ptr<Force1D>& f, L2_1D_EF v,
	VectorXd coeffs, Transform1D t)
{
	double err = 0.0;
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		double x = t.forwardTransform(x_hat);
		VectorXd v_vals = v.eval(x_hat)/t.jacobian(x_hat);
		err += pow(f->sol(x) - v_vals.dot(coeffs), 2)*t.jacobian(x_hat)*gl.getWeight(i);
	}
	return err;
}

class Integrator2D
{
	public:
		Integrator2D(json params);
		~Integrator2D() = default;
		MatrixXd laplace(H1_2D u, H1_2D v, Transform2D t);
		MatrixXd mass(HDiv_2D u, HDiv_2D v, Transform2D t);
		MatrixXd grad(L2_2D u, HDiv_2D v, Transform2D t);
		MatrixXd div(HDiv_2D u, L2_2D v, Transform2D t);
		MatrixXd mass(HCurl_2D u, HCurl_2D v, Transform2D t);
		MatrixXd grad(H1_2D u, HCurl_2D v, Transform2D t);
		MatrixXd div(HCurl_2D u, H1_2D v, Transform2D t);
		VectorXd force(const std::shared_ptr<Force2D>& f, H1_2D v, Transform2D t);
		VectorXd force(const std::shared_ptr<Force2D>& f, L2_2D v, Transform2D t);
		double error(const std::shared_ptr<Force2D>& f, H1_2D v,
			VectorXd coeffs, Transform2D t);
		double error(const std::shared_ptr<Force2D>& f, L2_2D v,
			VectorXd coeffs, Transform2D t);
	private:
		GaussLegendre gl;
};

Integrator2D::Integrator2D(json params) : gl((int)params["N"]) {}

inline MatrixXd Integrator2D::laplace(H1_2D u, H1_2D v, Transform2D t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int j = 0; j < gl.getN(); ++j)
	{
		Vector2d x_hat;
		x_hat[1] = gl.getNode(j);
		double w_j = gl.getWeight(j);
		for (int i = 0; i < gl.getN(); ++i)
		{
			x_hat[0] = gl.getNode(i);
			Matrix2d j_inv = t.jacobianInverse(x_hat);
			MatrixX2d u_vals = u.evalGrad(x_hat)*j_inv;
			MatrixX2d v_vals = v.evalGrad(x_hat)*j_inv;
			matrix += v_vals*u_vals.transpose()*t.jacobian(x_hat)*gl.getWeight(i)*w_j;
		}
	}
	return matrix;
}

inline MatrixXd Integrator2D::mass(HDiv_2D u, HDiv_2D v, Transform2D t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int j = 0; j < gl.getN(); ++j)
	{
		Vector2d x_hat;
		x_hat[1] = gl.getNode(j);
		double w_j = gl.getWeight(j);
		for (int i = 0; i < gl.getN(); ++i)
		{
			x_hat[0] = gl.getNode(i);
			Matrix2d j_t = t.jacobianMatrix(x_hat).transpose();
			MatrixX2d u_vals = u.eval(x_hat)*j_t;
			MatrixX2d v_vals = v.eval(x_hat)*j_t;
			matrix += v_vals*u_vals.transpose()/t.jacobian(x_hat)*gl.getWeight(i)*w_j;
		}
	}
	return matrix;
}

inline MatrixXd Integrator2D::grad(L2_2D u, HDiv_2D v, Transform2D t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int j = 0; j < gl.getN(); ++j)
	{
		Vector2d x_hat;
		x_hat[1] = gl.getNode(j);
		double w_j = gl.getWeight(j);
		for (int i = 0; i < gl.getN(); ++i)
		{
			x_hat[0] = gl.getNode(i);;
			VectorXd u_vals = u.eval(x_hat);
			VectorXd v_vals = v.evalDiv(x_hat);
			matrix += v_vals*u_vals.transpose()*gl.getWeight(i)*w_j;
		}
	}
	return matrix;
}

inline MatrixXd Integrator2D::div(HDiv_2D u, L2_2D v, Transform2D t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int j = 0; j < gl.getN(); ++j)
	{
		Vector2d x_hat;
		x_hat[1] = gl.getNode(j);
		double w_j = gl.getWeight(j);
		for (int i = 0; i < gl.getN(); ++i)
		{
			x_hat[0] = gl.getNode(i);;
			VectorXd u_vals = u.evalDiv(x_hat);
			VectorXd v_vals = v.eval(x_hat);
			matrix += v_vals*u_vals.transpose()*gl.getWeight(i)*w_j;
		}
	}
	return matrix;
}

inline MatrixXd Integrator2D::mass(HCurl_2D u, HCurl_2D v, Transform2D t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int j = 0; j < gl.getN(); ++j)
	{
		Vector2d x_hat;
		x_hat[1] = gl.getNode(j);
		double w_j = gl.getWeight(j);
		for (int i = 0; i < gl.getN(); ++i)
		{
			x_hat[0] = gl.getNode(i);
			Matrix2d j_inv = t.jacobianInverse(x_hat);
			MatrixX2d u_vals = u.eval(x_hat)*j_inv;
			MatrixX2d v_vals = v.eval(x_hat)*j_inv;
			matrix += v_vals*u_vals.transpose()*t.jacobian(x_hat)*gl.getWeight(i)*w_j;
		}
	}
	return matrix;
}

inline MatrixXd Integrator2D::grad(H1_2D u, HCurl_2D v, Transform2D t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int j = 0; j < gl.getN(); ++j)
	{
		Vector2d x_hat;
		x_hat[1] = gl.getNode(j);
		double w_j = gl.getWeight(j);
		for (int i = 0; i < gl.getN(); ++i)
		{
			x_hat[0] = gl.getNode(i);
			Matrix2d j_inv = t.jacobianInverse(x_hat);
			MatrixX2d u_vals = u.evalGrad(x_hat)*j_inv;
			MatrixX2d v_vals = v.eval(x_hat)*j_inv;
			matrix += v_vals*u_vals.transpose()*t.jacobian(x_hat)*gl.getWeight(i)*w_j;
		}
	}
	return matrix;
}

inline MatrixXd Integrator2D::div(HCurl_2D u, H1_2D v, Transform2D t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int j = 0; j < gl.getN(); ++j)
	{
		Vector2d x_hat;
		x_hat[1] = gl.getNode(j);
		double w_j = gl.getWeight(j);
		for (int i = 0; i < gl.getN(); ++i)
		{
			x_hat[0] = gl.getNode(i);
			Matrix2d j_inv = t.jacobianInverse(x_hat);
			MatrixX2d u_vals = u.eval(x_hat)*j_inv;
			MatrixX2d v_vals = v.evalGrad(x_hat)*j_inv;
			matrix += v_vals*u_vals.transpose()*t.jacobian(x_hat)*gl.getWeight(i)*w_j;
		}
	}
	return matrix;
}

inline VectorXd Integrator2D::force(const std::shared_ptr<Force2D>& f, H1_2D v,
	Transform2D t)
{
	VectorXd rhs(v.dofs());
	rhs.setZero();
	for (int j = 0; j < gl.getN(); ++j)
	{
		Vector2d x_hat;
		x_hat[1] = gl.getNode(j);
		double w_j = gl.getWeight(j);
		for (int i = 0; i < gl.getN(); ++i)
		{
			x_hat[0] = gl.getNode(i);
			Vector2d x = t.forwardTransform(x_hat);
			VectorXd v_vals = v.eval(x_hat);
			rhs += f->f(x)*v_vals*t.jacobian(x_hat)*gl.getWeight(i)*w_j;
		}
	}
	return rhs;
}

inline VectorXd Integrator2D::force(const std::shared_ptr<Force2D>& f, L2_2D v,
	Transform2D t)
{
	VectorXd rhs(v.dofs());
	rhs.setZero();
	for (int j = 0; j < gl.getN(); ++j)
	{
		Vector2d x_hat;
		x_hat[1] = gl.getNode(j);
		double w_j = gl.getWeight(j);
		for (int i = 0; i < gl.getN(); ++i)
		{
			x_hat[0] = gl.getNode(i);
			Vector2d x = t.forwardTransform(x_hat);
			VectorXd v_vals = v.eval(x_hat);
			rhs += f->f(x)*v_vals*t.jacobian(x_hat)*gl.getWeight(i)*w_j;
		}
	}
	return rhs;
}

inline double Integrator2D::error(const std::shared_ptr<Force2D>& f, H1_2D v,
	VectorXd coeffs, Transform2D t)
{
	double err = 0.0;
	for (int j = 0; j < gl.getN(); ++j)
	{
		Vector2d x_hat;
		x_hat[1] = gl.getNode(j);
		double w_j = gl.getWeight(j);
		for (int i = 0; i < gl.getN(); ++i)
		{
			x_hat[0] = gl.getNode(i);
			Vector2d x = t.forwardTransform(x_hat);
			VectorXd v_vals = v.eval(x_hat);
			err += pow(f->sol(x) - v_vals.dot(coeffs), 2)*t.jacobian(x_hat)*gl.getWeight(i)*w_j;
		}
	}
	return err;
}

inline double Integrator2D::error(const std::shared_ptr<Force2D>& f, L2_2D v,
	VectorXd coeffs, Transform2D t)
{
	double err = 0.0;
	for (int j = 0; j < gl.getN(); ++j)
	{
		Vector2d x_hat;
		x_hat[1] = gl.getNode(j);
		double w_j = gl.getWeight(j);
		for (int i = 0; i < gl.getN(); ++i)
		{
			x_hat[0] = gl.getNode(i);
			Vector2d x = t.forwardTransform(x_hat);
			VectorXd v_vals = v.eval(x_hat);
			err += pow(f->sol(x) - v_vals.dot(coeffs), 2)*t.jacobian(x_hat)*gl.getWeight(i)*w_j;
		}
	}
	return err;
}

#endif
