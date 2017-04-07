#include "json.hpp"
#include "ReferenceElements.hpp"
#include "ElementTransforms.hpp"

#ifndef _Integrators_hpp
#define _Integrators_hpp

using namespace nlohmann;
using namespace Eigen;

class Integrator1D
{
	public:
		Integrator1D(json params);
		~Integrator1D() = default;
		MatrixXd mass(H1_1D u, H1_1D v, Trasform1D_Linear t);
		MatrixXd mass(L2_1D u, L2_1D v, Trasform1D_Linear t);
		MatrixXd grad(H1_1D u, L2_1D v, Trasform1D_Linear t);
		MatrixXd grad(L2_1D u, H1_1D v, Trasform1D_Linear t);
		MatrixXd laplace(H1_1D u, H1_1D v, Trasform1D_Linear t);
		VectorXd force(const Force1D& f, H1_1D v, Trasform1D_Linear t);
		VectorXd force(const Force1D& f, L2_1D v, Trasform1D_Linear t);
	private:
		const GaussLegendre gl;
};

Integrator1D::Integrator1D(njson params) : gl((int)params["N"]) {}

inline MatrixXd Integrator1D::mass(H1_1D u, H1_1D v, Trasform1D_Linear t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd u_vals = u.eval(x_hat);
		VectorXd v_vals = v.eval(x_hat);
		matrix += v_vals*u_vals.traspose()*t.jacobian();
	}
	return matrix;
}

inline MatrixXd Integrator1D::mass(L2_1D u, L2_1D v, Trasform1D_Linear t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd u_vals = u.eval(x_hat);
		VectorXd v_vals = v.eval(x_hat);
		matrix += v_vals*u_vals.traspose()*t.jacobian();
	}
	return matrix;
}

inline MatrixXd Integrator1D::grad(H1_1D u, L2_1D v, Trasform1D_Linear t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd u_vals = u.evalD(x_hat);
		VectorXd v_vals = v.eval(x_hat);
		matrix += v_vals*u_vals.traspose();
	}
	return matrix;
}

inline MatrixXd Integrator1D::grad(L2_1D u, H1_1D v, Trasform1D_Linear t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd u_vals = u.eval(x_hat);
		VectorXd v_vals = v.evalD(x_hat);
		matrix += v_vals*u_vals.traspose();
	}
	return matrix;
}

inline MatrixXd Integrator1D::laplace(H1_1D u, H1_1D v, Trasform1D_Linear t)
{
	MatrixXd matrix(v.dofs(), u.dofs());
	matrix.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd u_vals = u.evalD(x_hat);
		VectorXd v_vals = v.evalD(x_hat);
		matrix += v_vals*u_vals.traspose()/t.jacobian();
	}
	return matrix;
}

inline VectorXd Integrator1D::force(const Force1D& f, H1_1D v, Trasform1D_Linear t)
{
	VectorXd rhs(v.dofs());
	rhs.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd v_vals = v.eval(x_hat);
		rhs += f(x_hat)*v_vals*t.jacobian();
	}
	return rhs;
}

inline VectorXd Integrator1D::force(const Force1D& f, L2_1D v, Trasform1D_Linear t)
{
	VectorXd rhs(v.dofs());
	rhs.setZero();
	for (int i = 0; i < gl.getN(); ++i)
	{
		double x_hat = gl.getNode(i);
		VectorXd v_vals = v.eval(x_hat);
		rhs += f(x_hat)*v_vals*t.jacobian();
	}
	return rhs;
}

#endif
