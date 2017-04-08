#include "json.hpp"
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include "Mesh1D.hpp"
#include "Integrators.hpp"
#include "ReferenceElements.hpp"

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
	H1_1D p(order);
	L2_1D u(order);
	std::shared_ptr<Force1D> force;

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

	std::cout << integrate.mass(u,u,mesh.getLinearTransform(0)) << std::endl;
	std::cout << std::endl;
	std::cout << integrate.grad(p,u,mesh.getLinearTransform(0)) << std::endl;
	output["x"] = mesh.nodes();

	std::ofstream outfile(argv[2]);
	outfile << std::setw(2) << output << std::endl;

	return 0;
}
