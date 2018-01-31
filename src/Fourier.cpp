#include "json.hpp"
#include <iostream>
#include "Integrators.hpp"
#include <complex>

using namespace Eigen;

int main(int argc, char const *argv[])
{
	int el_order = 1;
	if (argc > 2)
	{
		std::cerr << "Only one input value allowed." << std::endl;
		return -1;
	}
	el_order = std::atoi(argv[1]);

	json integrator_params = {
		{"N", 20}
	};
	Integrator2D integrator_2d(integrator_params);

	H1_2D h1_el(el_order);
	HCurl_2D hcurl_el(el_order);
	HDiv_2D hdiv_el(el_order);
	L2_2D l2_el(el_order);

	Transform2D transform_2d;

	MatrixXd grad_w_mass = integrator_2d.grad(h1_el, hcurl_el, transform_2d);
	MatrixXd mass_hcurl = integrator_2d.mass(hcurl_el, hcurl_el, transform_2d);

	std::cout << std::fixed << std::setprecision(3) << mass_hcurl << std::endl;

	MatrixXd mass_h1 = integrator_2d.mass(h1_el, h1_el, transform_2d);

	std::cout << std::fixed << std::setprecision(3) << mass_h1 << std::endl;

	MatrixXd grad = mass_hcurl.partialPivLu().solve(grad_w_mass);

	MatrixXd grad_3d(h1_el.dofs() + hcurl_el.dofs(), h1_el.dofs());
	grad_3d << grad, MatrixXd::Identity(h1_el.dofs(), h1_el.dofs());

	std::cout << std::fixed << std::setprecision(3) << grad_3d << std::endl;
	return 0;
}