#include "json.hpp"
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include "Eigen/Sparse"
#include "Mesh1D.hpp"
#include "Integrators.hpp"
#include "ReferenceElements.hpp"

using namespace Eigen;

int main(int argc, char const *argv[])
{
	if (argc != 3)
	{
		std::cerr << "Must provide the names of the input and output files." << std::endl;
		return 1;
	}

	json input;
	json output;

	std::ifstream infile(argv[1]);
	infile >> input;

	Mesh1D mesh(input["Mesh"]);
	Integrator1D integrate(input["Integrator"]);
	int h1_order = (int)input["H1_Order"];
	int l2_order = (int)input["L2_Order"];
	int N_el = (int)input["Mesh"]["N"] + 1;
	H1_1D h1(h1_order);
	L2_1D l2(l2_order);
	std::shared_ptr<Force1D> force;
	SparseMatrix<double> system;
	VectorXd rhs;
	// SimplicialLLT<SparseMatrix<double> > solver;
	SparseLU<SparseMatrix<double> > solver;

	std::string force_type = input["Function"];
	if (force_type == "Two")
	{
		force = std::make_shared<Two>();
	}
	else if (force_type == "ExpX3pX")
	{
		force = std::make_shared<ExpX3pX>();
	}
	else
	{
		std::cerr << "Unknown forcing function." << std::endl;
		return 2;
	}

	std::string formulation = input["Formulation"];
	if (formulation == "Standard")
	{
		// This is for the H1 standard system.
		system.resize(N_el*h1_order - 1, N_el*h1_order - 1);
		rhs.resize(N_el*h1_order - 1);
		system.reserve(VectorXi::Constant(N_el*h1_order - 1, 2*h1_order + 1));
		system.setZero();
		rhs.setZero();
		MatrixXd minimatrix = integrate.laplace(h1, h1, mesh.getTransform(0));
		VectorXd minirhs = integrate.force(force, h1, mesh.getTransform(0));
		for (int j = 1; j <= h1_order; ++j)
		{
			rhs[j-1] += minirhs[j];
			for (int i = 1; i <= h1_order; ++i)
			{
				system.coeffRef(i-1, j-1) += minimatrix(i, j);
			}
		}
		for (int n = 1; n < N_el - 1; ++n)
		{
			minimatrix = integrate.laplace(h1, h1, mesh.getTransform(n));
			minirhs = integrate.force(force, h1, mesh.getTransform(n));
			for (int j = 0; j <= h1_order; ++j)
			{
				rhs[n*h1_order+j-1] += minirhs[j];
				for (int i = 0; i <= h1_order; ++i)
				{
					system.coeffRef(n*h1_order+i-1, n*h1_order+j-1) += minimatrix(i, j);
				}
			}
		}
		minimatrix = integrate.laplace(h1, h1, mesh.getTransform(N_el - 1));
		minirhs = integrate.force(force, h1, mesh.getTransform(N_el - 1));
		for (int j = 0; j < h1_order; ++j)
		{
			rhs[(N_el-1)*h1_order+j-1] += minirhs[j];
			for (int i = 0; i < h1_order; ++i)
			{
				system.coeffRef((N_el-1)*h1_order+i-1,
					(N_el-1)*h1_order+j-1) += minimatrix(i, j);
			}
		}
	}
	else if(formulation == "Mixed")
	{
		/*
		 * (D) :
		 * 		u = -dp/dx
		 * 		du/dx = f
		 * 		p| = 0
		 * (V) :
		 * 	Look for u in H^1 and p in L^2 s.t.
		 * 		-<u,v> + <p,dv/dx> = 0	for all v in H^1
		 * 		<du/dx,q> = <f,q>		for all q in L^2
		 */
		system.resize(N_el*(h1_order + l2_order) + 1, N_el*(h1_order + l2_order) + 1);
		rhs.resize(N_el*(h1_order + l2_order) + 1);
		system.reserve(VectorXi::Constant(N_el*(h1_order + l2_order) + 1, 4*h1_order + 1));
		system.setZero();
		rhs.setZero();
		for (int n = 0; n < N_el; ++n)
		{
			MatrixXd mini_h1_mass = integrate.mass(h1, h1, mesh.getTransform(n));
			MatrixXd mini_div = integrate.grad(h1, l2, mesh.getTransform(n));
			for (int j = 0; j <= h1_order; ++j)
			{
				for (int i = 0; i <= h1_order; ++i)
				{
					system.coeffRef(n*h1_order+i, n*h1_order+j) -= mini_h1_mass(i, j);
				}
				for (int i = 0; i < l2_order; ++i)
				{
					system.coeffRef(N_el*h1_order+n*l2_order+i+1,
						n*h1_order+j) += mini_div(i, j);
				}
			}
			VectorXd minirhs = integrate.force(force, l2, mesh.getTransform(n));
			for (int i = 0; i < l2_order; ++i)
			{
				rhs[N_el*h1_order+n*l2_order+i+1] = minirhs[i];
			}
		}
		SparseMatrix<double> identity(N_el*(h1_order + l2_order) + 1,
			N_el*(h1_order + l2_order) + 1);
		identity.setIdentity();
		system = system.selfadjointView<Lower>()*identity;
	}
	else if (formulation == "Mimetic")
	{
		/*
		 * (D) :
		 * 		u = -dp/dx
		 * 		du/dx = f
		 * 		p| = 0
		 * (V) :
		 * 	Look for u in H^1 and p in L^2 s.t.
		 * 		<u,v> + <dp/dx,v> = 0	for all v in L^2
		 * 		<u,dq/dx> = -<f,q>		for all q in H^1
		 */
		system.resize(N_el*(h1_order + l2_order) - 1, N_el*(h1_order + l2_order) - 1);
		rhs.resize(N_el*(h1_order + l2_order) - 1);
		system.reserve(VectorXi::Constant(N_el*(h1_order + l2_order) - 1,
			l2_order + 2*h1_order));
		system.setZero();
		rhs.setZero();
		for (int n = 0; n < N_el; ++n)
		{
			MatrixXd mini_l2_mass = integrate.mass(l2, l2, mesh.getTransform(n));
			for (int j = 0; j < l2_order; ++j)
			{
				for (int i = 0; i < l2_order; ++i)
				{
					system.coeffRef(n*l2_order+i, n*l2_order+j) += mini_l2_mass(i, j);
				}
			}
		}
		MatrixXd mini_grad = integrate.grad(h1, l2, mesh.getTransform(0));
		VectorXd minirhs = integrate.force(force, h1, mesh.getTransform(0));
		for (int j = 1; j <= h1_order; ++j)
		{
			rhs[N_el*l2_order+j-1] -= minirhs[j];
			for (int i = 0; i < l2_order; ++i)
			{
				system.coeffRef(i, N_el*l2_order+j-1) += mini_grad(i, j);
			}
		}
		for (int n = 1; n < N_el - 1; ++n)
		{
			mini_grad = integrate.grad(h1, l2, mesh.getTransform(n));
			minirhs = integrate.force(force, h1, mesh.getTransform(n));
			for (int j = 0; j <= h1_order; ++j)
			{
				rhs[N_el*l2_order+n*h1_order+j-1] -= minirhs[j];
				for (int i = 0; i < l2_order; ++i)
				{
					system.coeffRef(n*l2_order+i,
						N_el*l2_order+n*h1_order+j-1) += mini_grad(i, j);
				}
			}
		}
		mini_grad = integrate.grad(h1, l2, mesh.getTransform(N_el - 1));
		minirhs = integrate.force(force, h1, mesh.getTransform(N_el - 1));
		for (int j = 0; j < h1_order; ++j)
		{
			rhs[N_el*l2_order+(N_el-1)*h1_order+j-1] -= minirhs[j];
			for (int i = 0; i < l2_order; ++i)
			{
				system.coeffRef((N_el-1)*l2_order+i,
					N_el*l2_order+(N_el-1)*h1_order+j-1) += mini_grad(i, j);
			}
		}
		SparseMatrix<double> identity(N_el*(h1_order + l2_order) - 1,
			N_el*(h1_order + l2_order) - 1);
		identity.setIdentity();
		system = system.selfadjointView<Upper>()*identity;
	}
	else
	{
		std::cerr << "Unknown formulation." << std::endl;
		return 3;
	}


	system.makeCompressed();

	std::cout << system << std::endl;
	std::cout << rhs << std::endl;
	std::cout << std::endl;

	solver.analyzePattern(system);
	solver.factorize(system);
	// std::cout << solver.determinant() << std::endl;
	// std::cout << std::endl;
	VectorXd x = solver.solve(rhs);

	std::cout << x << std::endl;

	double error = 0.0;
	if (formulation == "Standard")
	{
		VectorXd segment(h1_order + 1);
		segment[0] = 0.0;
		segment.tail(h1_order) = x.segment(0, h1_order);
		error += integrate.error(force, h1, segment, mesh.getTransform(0));
		for (int n = 1; n < N_el - 1; ++n)
		{
			error += integrate.error(force, h1,
				x.segment(n*h1_order - 1, h1_order + 1), mesh.getTransform(n));
		}
		segment.head(h1_order) = x.segment((N_el-1)*h1_order-1, h1_order);
		segment[h1_order] = 0.0;
		error += integrate.error(force, h1, segment, mesh.getTransform(N_el-1));
	}
	else if(formulation == "Mixed")
	{
		for (int n = 0; n < N_el; ++n)
		{
			error += integrate.error(force, l2,
				x.segment(N_el*h1_order+n*l2_order+1, l2_order), mesh.getTransform(n));
		}
	}
	else if (formulation == "Mimetic")
	{
		VectorXd segment(h1_order + 1);
		segment[0] = 0.0;
		segment.tail(h1_order) = x.segment(N_el*l2_order, h1_order);
		error += integrate.error(force, h1, segment, mesh.getTransform(0));
		for (int n = 1; n < N_el - 1; ++n)
		{
			error += integrate.error(force, h1,
				x.segment(N_el*l2_order+n*h1_order-1, h1_order + 1), mesh.getTransform(n));
		}
		segment.head(h1_order) = x.segment(N_el*l2_order+(N_el-1)*h1_order-1, h1_order);
		segment[h1_order] = 0.0;
		error += integrate.error(force, h1, segment, mesh.getTransform(N_el-1));
	}

	error = sqrt(error);

	std::cout << std::endl;
	std::cout << error << std::endl;

	// // std::vector<double> v(N_el*order + 1);
	// std::vector<double> v(N_el*order);
	// // v[0] = 0.0;
	// v[0] = x[N_el*order-2];
	// for (int i = 1; i < N_el*order; ++i)
	// {
	// 	// v[i] = x[i-1];
	// 	v[i] = x[N_el*order+i-2];
	// }
	// // v[N_el*order] = 0.0;

	// std::cout << integrate.mass(l2, l2, mesh.getTransform(0)) << std::endl;
	// std::cout << std::endl;
	// std::cout << integrate.grad(l2, h1, mesh.getTransform(0)) << std::endl;
	// std::cout << std::endl;
	// std::cout << integrate.force(force, h1, mesh.getTransform(0)) << std::endl;
	// std::cout << std::endl;
	// std::cout << integrate.mass(h1, h1, mesh.getTransform(0)) << std::endl;
	// std::cout << std::endl;
	// std::cout << integrate.grad(h1,l2, mesh.getTransform(0)) << std::endl;
	// std::cout << std::endl;
	// std::cout << integrate.force(force, l2, mesh.getTransform(0)) << std::endl;

	output["x"] = mesh.nodes();
	output["error"] = error;
	// output["u"] = v;

	std::ofstream outfile(argv[2]);
	outfile << std::setw(2) << output << std::endl;

	return 0;
}
