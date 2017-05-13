#include "json.hpp"
#include <iostream>
#include <fstream>
#include "Integrators.hpp"
#include "ForcingFunctions.hpp"
#include "Eigen/Eigenvalues"
#include "Mesh2D.hpp"

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

	Mesh2D mesh(input["Mesh"]);
	Integrator2D integrate(input["Integrator"]);
	std::shared_ptr<Force2D> force;
	SparseMatrix<double> system;
	VectorXd rhs;
	SparseLU<SparseMatrix<double> > solver;

	std::string force_type = input["Function"];
	if (force_type == "SinXSinY")
	{
		force = std::make_shared<SinXSinY>();
	}
	else
	{
		std::cerr << "Unknown forcing function." << std::endl;
		return 2;
	}

	std::string formulation = input["Formulation"];
	if (formulation == "Standard")
	{
		system = mesh.assembleStandard(integrate);
		rhs = mesh.rhsStandard(integrate, force);
	}
	else if(formulation == "Mixed")
	{
		system = mesh.assembleMixed(integrate);
		rhs = mesh.rhsMixed(integrate, force);
	}
	else if(formulation == "Dual-Mixed")
	{
		system = mesh.assembleDualMixed(integrate);
		rhs = mesh.rhsDualMixed(integrate, force);
	}
	else
	{
		std::cerr << "Unknown formulation." << std::endl;
		return 3;
	}

	solver.analyzePattern(system);
	solver.factorize(system);
	VectorXd x = solver.solve(rhs);

	// std::cout << x << std::endl;

	double error;
	if (formulation == "Standard")
	{
		error = mesh.errorStandard(integrate, force, x);
	}
	else if(formulation == "Mixed")
	{
		error = mesh.errorMixed(integrate, force, x);
	}
	else if(formulation == "Dual-Mixed")
	{
		error = mesh.errorDualMixed(integrate, force, x);
	}

	error = sqrt(error);
	std::cout << error << std::endl;

	return 0;
}