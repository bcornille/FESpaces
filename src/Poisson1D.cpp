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
	int order = (int)input["Order"];
	int N_el = (int)input["Mesh"]["N"] + 1;
	H1_1D h1(order);
	L2_1D l2(order);
	std::shared_ptr<Force1D> force;
	SparseMatrix<double> system;
	VectorXd rhs;
	SimplicialLLT<SparseMatrix<double> > solver;

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
		system.resize(N_el*order - 1, N_el*order - 1);
		rhs.resize(N_el*order - 1);
		system.reserve(VectorXi::Constant(N_el*order - 1, 2*order + 1));
		system.setZero();
		rhs.setZero();
		MatrixXd minimatrix = integrate.laplace(h1, h1, mesh.getLinearTransform(0));
		VectorXd minirhs = integrate.force(force, h1, mesh.getLinearTransform(0));
		for (int j = 1; j <= order; ++j)
		{
			rhs[j-1] += minirhs[j];
			for (int i = 1; i <= order; ++i)
			{
				system.coeffRef(i-1, j-1) += minimatrix(i, j);
			}
		}
		for (int n = 1; n < N_el - 1; ++n)
		{
			minimatrix = integrate.laplace(h1, h1, mesh.getLinearTransform(n));
			minirhs = integrate.force(force, h1, mesh.getLinearTransform(n));
			for (int j = 0; j <= order; ++j)
			{
				rhs[n*order+j-1] += minirhs[j];
				for (int i = 0; i <= order; ++i)
				{
					system.coeffRef(n*order+i-1, n*order+j-1) += minimatrix(i, j);
				}
			}
		}
		minimatrix = integrate.laplace(h1, h1, mesh.getLinearTransform(N_el - 1));
		minirhs = integrate.force(force, h1, mesh.getLinearTransform(N_el - 1));
		for (int j = 0; j < order; ++j)
		{
			rhs[(N_el-1)*order+j-1] += minirhs[j];
			for (int i = 0; i < order; ++i)
			{
				system.coeffRef((N_el-1)*order+i-1, (N_el-1)*order+j-1) += minimatrix(i, j);
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
		system.resize(2*N_el*order - 1, 2*N_el*order - 1);
		rhs.resize(2*N_el*order - 1);
		system.reserve(VectorXi::Constant(2*N_el*order - 1, 3*order + 1));
		system.setZero();
		rhs.setZero();
		MatrixXd mini_h1_mass = integrate.mass(h1, h1, mesh.getLinearTransform(0));
		MatrixXd mini_div = integrate.grad(h1, l2, mesh.getLinearTransform(0));
		VectorXd minirhs = integrate.force(force, l2, mesh.getLinearTransform(0));
		for (int j = 1; j <= order; ++j)
		{
			for (int i = 1; i <= order; ++i)
			{
				system.coeffRef(i-1, j-1) -= mini_h1_mass(i, j);
			}
			for (int i = 0; i < order; ++i)
			{
				system.coeffRef(N_el*order+i-1, j-1) += mini_div(i, j);
			}
		}
		for (int i = 0; i < order; ++i)
		{
			rhs[N_el*order+i-1] = minirhs[i];
		}
		for (int n = 1; n < N_el - 1; ++n)
		{
			mini_h1_mass = integrate.mass(h1, h1, mesh.getLinearTransform(n));
			mini_div = integrate.grad(h1, l2, mesh.getLinearTransform(n));
			for (int j = 0; j <= order; ++j)
			{
				for (int i = 0; i <= order; ++i)
				{
					system.coeffRef(n*order+i-1, n*order+j-1) -= mini_h1_mass(i, j);
				}
				for (int i = 0; i < order; ++i)
				{
					system.coeffRef((N_el+n)*order+i-1, n*order+j-1) += mini_div(i, j);
				}
			}
			for (int i = 0; i < order; ++i)
			{
				rhs[(N_el+n)*order+i-1] = minirhs[i];
			}
		}
		mini_h1_mass = integrate.mass(h1, h1, mesh.getLinearTransform(N_el-1));
		mini_div = integrate.grad(h1, l2, mesh.getLinearTransform(N_el-1));
		for (int j = 0; j < order; ++j)
		{
			for (int i = 0; i < order; ++i)
			{
				system.coeffRef((N_el-1)*order+i-1, (N_el-1)*order+j-1) -= mini_h1_mass(i, j);
			}
			for (int i = 0; i < order; ++i)
			{
				system.coeffRef((2*N_el-1)*order+i-1, (N_el-1)*order+j-1) += mini_div(i, j);
			}
		}
		for (int i = 0; i < order; ++i)
		{
			rhs[(2*N_el-1)*order+i-1] = minirhs[i];
		}
		SparseMatrix<double> identity(2*N_el*order - 1, 2*N_el*order - 1);
		identity.setIdentity();
		system = system.selfadjointView<Lower>()*identity;
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
	std::cout << solver.determinant() << std::endl;
	std::cout << std::endl;
	VectorXd x = solver.solve(rhs);

	std::cout << x << std::endl;

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

	// std::cout << integrate.mass(l2, l2, mesh.getLinearTransform(0)) << std::endl;
	// std::cout << std::endl;
	// std::cout << integrate.grad(l2, h1, mesh.getLinearTransform(0)) << std::endl;
	// std::cout << std::endl;
	// std::cout << integrate.force(force, h1, mesh.getLinearTransform(0)) << std::endl;
	// std::cout << std::endl;
	// std::cout << integrate.mass(h1, h1, mesh.getLinearTransform(0)) << std::endl;
	// std::cout << std::endl;
	// std::cout << integrate.grad(h1,l2, mesh.getLinearTransform(0)) << std::endl;
	// std::cout << std::endl;
	// std::cout << integrate.force(force, l2, mesh.getLinearTransform(0)) << std::endl;

	output["x"] = mesh.nodes();
	// output["u"] = v;

	std::ofstream outfile(argv[2]);
	outfile << std::setw(2) << output << std::endl;

	return 0;
}
