#include "json.hpp"
#include "ReferenceElements.hpp"
#include "ElementTransforms.hpp"
#include "ForcingFunctions.hpp"

#ifndef _Integrators_hpp
#define _Integrators_hpp

using namespace nlohmann;
using namespace Eigen;

class Integrator1D
{
	public:
		Integrator1D(json params);
		~Integrator1D() = default;
		MatrixXd mass(H1_1D u, H1_1D v, Transform1D_Linear t);
		MatrixXd mass(L2_1D u, L2_1D v, Transform1D_Linear t);
		MatrixXd mass(L2_1D_EF u, L2_1D_EF v, Transform1D_Linear t);
		MatrixXd grad(H1_1D u, L2_1D v, Transform1D_Linear t);
		MatrixXd grad(L2_1D u, H1_1D v, Transform1D_Linear t);
		MatrixXd grad(H1_1D u, L2_1D_EF v, Transform1D_Linear t);
		MatrixXd grad(L2_1D_EF u, H1_1D v, Transform1D_Linear t);
		MatrixXd laplace(H1_1D u, H1_1D v, Transform1D_Linear t);
		VectorXd force(const std::shared_ptr<Force1D>& f, H1_1D v, Transform1D_Linear t);
		VectorXd force(const std::shared_ptr<Force1D>& f, L2_1D v, Transform1D_Linear t);
		VectorXd force(const std::shared_ptr<Force1D>& f, L2_1D_EF v, Transform1D_Linear t);
		double error(const std::shared_ptr<Force1D>& f, H1_1D v, VectorXd coeffs, Transform1D_Linear t);
		double error(const std::shared_ptr<Force1D>& f, L2_1D v, VectorXd coeffs, Transform1D_Linear t);
		double error(const std::shared_ptr<Force1D>& f, L2_1D_EF v, VectorXd coeffs, Transform1D_Linear t);
	private:
		GaussLegendre gl;
};

Integrator1D::Integrator1D(json params) : gl((int)params["N"]) {}

inline MatrixXd Integrator1D::mass(H1_1D u, H1_1D v, Transform1D_Linear t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd u_vals = u.eval(x_hat);
		VectorXd v_vals = v.eval(x_hat);
		matrix += v_vals*u_vals.transpose()*t.jacobian()*gl.getWeight(i);
	}
	return matrix;
}

inline MatrixXd Integrator1D::mass(L2_1D u, L2_1D v, Transform1D_Linear t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd u_vals = u.eval(x_hat);
		VectorXd v_vals = v.eval(x_hat);
		matrix += v_vals*u_vals.transpose()*t.jacobian()*gl.getWeight(i);
	}
	return matrix;
}

inline MatrixXd Integrator1D::mass(L2_1D_EF u, L2_1D_EF v, Transform1D_Linear t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd u_vals = u.eval(x_hat)/t.jacobian();
		VectorXd v_vals = v.eval(x_hat)/t.jacobian();
		matrix += v_vals*u_vals.transpose()*t.jacobian()*gl.getWeight(i);
	}
	return matrix;
}

inline MatrixXd Integrator1D::grad(H1_1D u, L2_1D v, Transform1D_Linear t)
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

inline MatrixXd Integrator1D::grad(L2_1D u, H1_1D v, Transform1D_Linear t)
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

inline MatrixXd Integrator1D::grad(H1_1D u, L2_1D_EF v, Transform1D_Linear t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd u_vals = u.evalD(x_hat);
		VectorXd v_vals = v.eval(x_hat)/t.jacobian();
		matrix += v_vals*u_vals.transpose()*gl.getWeight(i);
	}
	return matrix;
}

inline MatrixXd Integrator1D::grad(L2_1D_EF u, H1_1D v, Transform1D_Linear t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd u_vals = u.eval(x_hat)/t.jacobian();
		VectorXd v_vals = v.evalD(x_hat);
		matrix += v_vals*u_vals.transpose()*gl.getWeight(i);
	}
	return matrix;
}

inline MatrixXd Integrator1D::laplace(H1_1D u, H1_1D v, Transform1D_Linear t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd u_vals = u.evalD(x_hat);
		VectorXd v_vals = v.evalD(x_hat);
		matrix += v_vals*u_vals.transpose()/t.jacobian()*gl.getWeight(i);
	}
	return matrix;
}

inline VectorXd Integrator1D::force(const std::shared_ptr<Force1D>& f, H1_1D v,
	Transform1D_Linear t)
{
	VectorXd rhs(v.dofs());
	rhs.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		double x = t.forwardTransform(x_hat);
		VectorXd v_vals = v.eval(x_hat);
		rhs += f->f(x)*v_vals*t.jacobian()*gl.getWeight(i);
	}
	return rhs;
}

inline VectorXd Integrator1D::force(const std::shared_ptr<Force1D>& f, L2_1D v,
	Transform1D_Linear t)
{
	VectorXd rhs(v.dofs());
	rhs.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		double x = t.forwardTransform(x_hat);
		VectorXd v_vals = v.eval(x_hat);
		rhs += f->f(x)*v_vals*t.jacobian()*gl.getWeight(i);
	}
	return rhs;
}

inline VectorXd Integrator1D::force(const std::shared_ptr<Force1D>& f, L2_1D_EF v,
	Transform1D_Linear t)
{
	VectorXd rhs(v.dofs());
	rhs.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		double x = t.forwardTransform(x_hat);
		VectorXd v_vals = v.eval(x_hat)/t.jacobian();
		rhs += f->f(x)*v_vals*t.jacobian()*gl.getWeight(i);
	}
	return rhs;
}

inline double Integrator1D::error(const std::shared_ptr<Force1D>& f, H1_1D v,
	VectorXd coeffs, Transform1D_Linear t)
{
	double err = 0.0;
	int N = gl.getN();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		double x = t.forwardTransform(x_hat);
		VectorXd v_vals = v.eval(x_hat);
		err += pow(f->sol(x) - v_vals.dot(coeffs), 2)*t.jacobian()*gl.getWeight(i);
	}
	return err;
}

inline double Integrator1D::error(const std::shared_ptr<Force1D>& f, L2_1D v,
	VectorXd coeffs, Transform1D_Linear t)
{
	double err = 0.0;
	int N = gl.getN();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		double x = t.forwardTransform(x_hat);
		VectorXd v_vals = v.eval(x_hat);
		err += pow(f->sol(x) - v_vals.dot(coeffs), 2)*t.jacobian()*gl.getWeight(i);
	}
	return err;
}

inline double Integrator1D::error(const std::shared_ptr<Force1D>& f, L2_1D_EF v,
	VectorXd coeffs, Transform1D_Linear t)
{
	double err = 0.0;
	int N = gl.getN();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		double x = t.forwardTransform(x_hat);
		VectorXd v_vals = v.eval(x_hat)/t.jacobian();
		err += pow(f->sol(x) - v_vals.dot(coeffs), 2)*t.jacobian()*gl.getWeight(i);
	}
	return err;
}

#endif
