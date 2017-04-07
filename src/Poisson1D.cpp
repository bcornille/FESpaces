#include "json.hpp"
#include <iostream>
#include <fstream>
#include <vector>
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

	output["x"] = mesh.nodes();

	std::ofstream outfile(argv[2]);
	outfile << std::setw(2) << output << std::endl;

	return 0;
}
