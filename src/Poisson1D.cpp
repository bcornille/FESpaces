#include "json.hpp"
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include "Mesh1D.hpp"
#include "Integrators.hpp"

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
	std::unique_ptr<Force1D> force;

	std::string force_type = input["Function"];
	if (force_type == "Two")
	{
		force = std::make_unique<Force1D>(new Two);
	}
	else if (force_type == "ExpX1mX")
	{
		force = std::make_unique<Force1D>(new ExpX1mX);
	}
	else
	{
		std::cerr << "Unknown forcing function." << std::endl;
		return 2;
	}

	output["x"] = mesh.nodes();

	std::ofstream outfile(argv[2]);
	outfile << std::setw(2) << output << std::endl;

	return 0;
}
