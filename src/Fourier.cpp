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
	if (argc == 2) el_order = std::atoi(argv[1]);

	json integrator_params = {
		{"N", 20}
	};
	Integrator2D integrator_2d(integrator_params);

	H1_2D h1_el(el_order);
	HCurl_2D hcurl_el(el_order);
	HDiv_2D hdiv_el(el_order);
	L2_2D l2_el(el_order);

	Transform2D transform_2d(Vector2d::Constant(-1.0), Vector2d::Constant(1.0));

	MatrixXd grad_w_mass = integrator_2d.grad(h1_el, hcurl_el, transform_2d);
	MatrixXd mass_hcurl = integrator_2d.mass(hcurl_el, hcurl_el, transform_2d);

	// std::cout << std::fixed << std::setprecision(3) << mass_hcurl << std::endl;

	MatrixXcd mass_h1 = integrator_2d.mass(h1_el, h1_el, transform_2d);

	// std::cout << std::fixed << std::setprecision(3) << mass_h1 << std::endl;

	MatrixXd grad = mass_hcurl.partialPivLu().solve(grad_w_mass);

	MatrixXcd grad_3d(h1_el.dofs() + hcurl_el.dofs(), h1_el.dofs());
	grad_3d << grad, 1.0i*MatrixXcd::Identity(h1_el.dofs(), h1_el.dofs());

	std::cout << std::fixed << std::setprecision(3) << grad_3d << std::endl;
	std::cout << std::endl;

	MatrixXcd grad_3d_w_mass(h1_el.dofs() + hcurl_el.dofs(), h1_el.dofs());
	grad_3d_w_mass << grad_w_mass, 1.0i*mass_h1;

	std::cout << std::fixed << std::setprecision(3) << grad_3d_w_mass << std::endl;
	std::cout << std::endl;

	MatrixXd curl_w_mass = integrator_2d.curl(hcurl_el, l2_el, transform_2d);
	MatrixXd mass_l2 = integrator_2d.mass(l2_el, l2_el, transform_2d);

	// std::cout << std::fixed << std::setprecision(3) << curl_w_mass << std::endl;
	// std::cout << std::endl;
	// std::cout << std::fixed << std::setprecision(3) << mass_l2 << std::endl;
	// std::cout << std::endl;

	MatrixXd curl = mass_l2.partialPivLu().solve(curl_w_mass);

	// std::cout << std::fixed << std::setprecision(3) << curl << std::endl;
	// std::cout << std::endl;

	MatrixXd rot_w_mass = integrator_2d.curl(h1_el, hdiv_el, transform_2d);
	MatrixXd mass_hdiv = integrator_2d.mass(hdiv_el, hdiv_el, transform_2d);

	// std::cout << std::fixed << std::setprecision(3) << rot_w_mass << std::endl;
	// std::cout << std::endl;
	// std::cout << std::fixed << std::setprecision(3) << mass_hdiv << std::endl;

	MatrixXd rot = mass_hdiv.partialPivLu().solve(rot_w_mass);

	// std::cout << std::fixed << std::setprecision(3) << rot << std::endl;
	// std::cout << std::endl;

	int n_dofs = hdiv_el.dofs();
	MatrixXcd ddphi(n_dofs, n_dofs);
	ddphi << 1.0i*MatrixXcd::Identity(n_dofs/2, n_dofs/2), MatrixXcd::Zero(n_dofs/2, n_dofs/2),
			MatrixXcd::Zero(n_dofs/2, n_dofs/2), -1.0i*MatrixXcd::Identity(n_dofs/2, n_dofs/2);

	// std::cout << std::fixed << std::setprecision(3) << ddphi << std::endl;

	MatrixXcd curl_3d(n_dofs + l2_el.dofs(), n_dofs + h1_el.dofs());
	curl_3d << ddphi, rot,
				curl, MatrixXcd::Zero(l2_el.dofs(), h1_el.dofs());

	std::cout << std::fixed << std::setprecision(3) << curl_3d << std::endl;
	std::cout << std::endl;

	MatrixXcd curl_3d_w_mass(n_dofs + l2_el.dofs(), n_dofs + h1_el.dofs());
	curl_3d_w_mass << mass_hdiv*ddphi, rot_w_mass,
						curl_w_mass, MatrixXcd::Zero(l2_el.dofs(), h1_el.dofs());

	std::cout << std::fixed << std::setprecision(3) << curl_3d_w_mass << std::endl;
	std::cout << std::endl;

	MatrixXd div_w_mass = integrator_2d.div(hdiv_el, l2_el, transform_2d);
	MatrixXd divg = mass_l2.partialPivLu().solve(div_w_mass);

	MatrixXcd div_3d(l2_el.dofs(), n_dofs + l2_el.dofs());
	div_3d << divg, 1.0i*MatrixXcd::Identity(l2_el.dofs(), l2_el.dofs());

	std::cout << std::fixed << std::setprecision(3) << div_3d << std::endl;
	std::cout << std::endl;

	MatrixXcd div_3d_w_mass(l2_el.dofs(), n_dofs + l2_el.dofs());
	MatrixXcd mass_l2_cd = mass_l2;
	div_3d_w_mass << div_w_mass, 1.0i*mass_l2_cd;

	std::cout << std::fixed << std::setprecision(3) << div_3d_w_mass << std::endl;
	std::cout << std::endl;

	return 0;
}