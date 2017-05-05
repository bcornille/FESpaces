#include "json.hpp"
#include <iostream>
#include <fstream>
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

	Integrator2D integrate(input["Integrator"]);
	int order = (int)input["Order"];
	H1_2D h1(order);
	HDiv_2D hdiv(order);
	L2_2D l2(order);
	HCurl_2D hcurl(order);
	Transform2D t(Eigen::Vector2d::Constant(0.0), Eigen::Vector2d::Constant(1.0));

	// std::cout << integrate.laplace(h1, h1, t) << std::endl;
	// std::cout << std::endl;
	// std::cout << integrate.mass(hdiv, hdiv, t) << std::endl;
	// std::cout << std::endl;
	// std::cout << integrate.grad(l2, hdiv, t) << std::endl;
	// std::cout << std::endl;
	// std::cout << integrate.div(hdiv, l2, t) << std::endl;
	// td::cout << std::endl;
	// std::cout << integrate.mass(hcurl, hcurl, t) << std::endl;
	// std::cout << std::endl;
	// std::cout << integrate.grad(h1, hcurl, t) << std::endl;
	// std::cout << std::endl;
	// std::cout << integrate.mass(hcurl, hcurl, t).partialPivLu().solve(integrate.grad(h1, hcurl, t)) << std::endl;

	return 0;
}