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
	H1_1D p(order);
	L2_1D_EF u(order);
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

	// This is for the H1 standard system.
	system.resize(N_el*order - 1, N_el*order - 1);
	rhs.resize(N_el*order - 1);
	system.reserve(VectorXi::Constant(N_el*order - 1, 2*order + 1));
	system.setZero();
	rhs.setZero();
	MatrixXd minimatrix = integrate.laplace(p, p, mesh.getLinearTransform(0));
	VectorXd minirhs = integrate.force(force, p, mesh.getLinearTransform(0));
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
		minimatrix = integrate.laplace(p, p, mesh.getLinearTransform(n));
		minirhs = integrate.force(force, p, mesh.getLinearTransform(n));
		for (int j = 0; j <= order; ++j)
		{
			rhs[n*order+j-1] += minirhs[j];
			for (int i = 0; i <= order; ++i)
			{
				system.coeffRef(n*order+i-1, n*order+j-1) += minimatrix(i, j);
			}
		}
	}
	minimatrix = integrate.laplace(p, p, mesh.getLinearTransform(N_el - 1));
	minirhs = integrate.force(force, p, mesh.getLinearTransform(N_el - 1));
	for (int j = 0; j < order; ++j)
	{
		rhs[(N_el-1)*order+j-1] += minirhs[j];
		for (int i = 0; i < order; ++i)
		{
			system.coeffRef((N_el-1)*order+i-1, (N_el-1)*order+j-1) += minimatrix(i, j);
		}
	}
	system.makeCompressed();
	solver.analyzePattern(system);
	solver.factorize(system);
	VectorXd x = solver.solve(rhs);

	std::vector<double> v(N_el*order + 1);
	v[0] = 0.0;
	for (int i = 1; i < N_el*order; ++i)
	{
		v[i] = x[i-1];
	}
	v[N_el*order] = 0.0;

	// std::cout << integrate.mass(u, u, mesh.getLinearTransform(0)) << std::endl;
	// std::cout << std::endl;
	// std::cout << integrate.grad(u, p, mesh.getLinearTransform(0)) << std::endl;
	// std::cout << std::endl;
	std::cout << integrate.force(force, p, mesh.getLinearTransform(0)) << std::endl;
	// std::cout << std::endl;
	// std::cout << integrate.mass(p, p, mesh.getLinearTransform(0)) << std::endl;
	// std::cout << std::endl;
	// std::cout << integrate.grad(p, u, mesh.getLinearTransform(0)) << std::endl;
	// std::cout << std::endl;
	// std::cout << integrate.force(force, u, mesh.getLinearTransform(0)) << std::endl;

	output["x"] = mesh.nodes();
	output["u"] = v;

	std::ofstream outfile(argv[2]);
	outfile << std::setw(2) << output << std::endl;

	return 0;
}
