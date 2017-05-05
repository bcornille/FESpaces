#include "json.hpp"
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include "Eigen/Sparse"
#include "Mesh1D.hpp"
#include "Integrators.hpp"
#include "ReferenceElements.hpp"
#include <cstdlib>

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
	int size = 0;
	H1_1D h1(order);
	L2_1D l2(order);
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
		size = N_el*order - 1;
		system.resize(size, size);
		rhs.resize(size);
		system.reserve(VectorXi::Constant(N_el*order - 1, 2*order + 1));
		system.setZero();
		rhs.setZero();
		MatrixXd minimatrix = integrate.laplace(h1, h1, mesh.getTransform(0));
		VectorXd minirhs = integrate.force(force, h1, mesh.getTransform(0));
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
			minimatrix = integrate.laplace(h1, h1, mesh.getTransform(n));
			minirhs = integrate.force(force, h1, mesh.getTransform(n));
			for (int j = 0; j <= order; ++j)
			{
				rhs[n*order+j-1] += minirhs[j];
				for (int i = 0; i <= order; ++i)
				{
					system.coeffRef(n*order+i-1, n*order+j-1) += minimatrix(i, j);
				}
			}
		}
		minimatrix = integrate.laplace(h1, h1, mesh.getTransform(N_el - 1));
		minirhs = integrate.force(force, h1, mesh.getTransform(N_el - 1));
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
		size = 2*N_el*order + 1;
		system.resize(size, size);
		rhs.resize(size);
		system.reserve(VectorXi::Constant(2*N_el*order + 1, 4*order + 1));
		system.setZero();
		rhs.setZero();
		for (int n = 0; n < N_el; ++n)
		{
			MatrixXd mini_h1_mass = integrate.mass(h1, h1, mesh.getTransform(n));
			MatrixXd mini_div = integrate.grad(h1, l2, mesh.getTransform(n));
			for (int j = 0; j <= order; ++j)
			{
				for (int i = 0; i <= order; ++i)
				{
					system.coeffRef(n*order+i, n*order+j) -= mini_h1_mass(i, j);
				}
				for (int i = 0; i < order; ++i)
				{
					system.coeffRef((N_el+n)*order+i+1, n*order+j) += mini_div(i, j);
				}
			}
			VectorXd minirhs = integrate.force(force, l2, mesh.getTransform(n));
			for (int i = 0; i < order; ++i)
			{
				rhs[(N_el+n)*order+i+1] = minirhs[i];
			}
		}
		SparseMatrix<double> identity(2*N_el*order + 1, 2*N_el*order + 1);
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
		size = 2*N_el*order - 1;
		system.resize(size, size);
		rhs.resize(size);
		system.reserve(VectorXi::Constant(2*N_el*order - 1, 3*order));
		system.setZero();
		rhs.setZero();
		for (int n = 0; n < N_el; ++n)
		{
			MatrixXd mini_l2_mass = integrate.mass(l2, l2, mesh.getTransform(n));
			for (int j = 0; j < order; ++j)
			{
				for (int i = 0; i < order; ++i)
				{
					system.coeffRef(n*order+i, n*order+j) += mini_l2_mass(i, j);
				}
			}
		}
		MatrixXd mini_grad = integrate.grad(h1, l2, mesh.getTransform(0));
		VectorXd minirhs = integrate.force(force, h1, mesh.getTransform(0));
		for (int j = 1; j <= order; ++j)
		{
			rhs[N_el*order+j-1] -= minirhs[j];
			for (int i = 0; i < order; ++i)
			{
				system.coeffRef(i, N_el*order+j-1) += mini_grad(i, j);
			}
		}
		for (int n = 1; n < N_el - 1; ++n)
		{
			mini_grad = integrate.grad(h1, l2, mesh.getTransform(n));
			minirhs = integrate.force(force, h1, mesh.getTransform(n));
			for (int j = 0; j <= order; ++j)
			{
				rhs[(N_el+n)*order+j-1] -= minirhs[j];
				for (int i = 0; i < order; ++i)
				{
					system.coeffRef(n*order+i, (N_el+n)*order+j-1) += mini_grad(i, j);
				}
			}
		}
		mini_grad = integrate.grad(h1, l2, mesh.getTransform(N_el - 1));
		minirhs = integrate.force(force, h1, mesh.getTransform(N_el - 1));
		for (int j = 0; j < order; ++j)
		{
			rhs[(2*N_el-1)*order+j-1] -= minirhs[j];
			for (int i = 0; i < order; ++i)
			{
				system.coeffRef((N_el-1)*order+i, (2*N_el-1)*order+j-1) += mini_grad(i, j);
			}
		}
		SparseMatrix<double> identity(2*N_el*order - 1, 2*N_el*order - 1);
		identity.setIdentity();
		system = system.selfadjointView<Upper>()*identity;
	}
	else
	{
		std::cerr << "Unknown formulation." << std::endl;
		return 3;
	}


	system.makeCompressed();

	// std::cout << system << std::endl;
	// std::cout << rhs << std::endl;
	// std::cout << std::endl;

	solver.analyzePattern(system);
	solver.factorize(system);
	// std::cout << solver.determinant() << std::endl;
	// std::cout << std::endl;
	VectorXd x = solver.solve(rhs);

	// std::cout << x << std::endl;

	double error = 0.0;
	if (formulation == "Standard")
	{
		VectorXd segment(order + 1);
		segment[0] = 0.0;
		segment.tail(order) = x.segment(0, order);
		error += integrate.error(force, h1, segment, mesh.getTransform(0));
		for (int n = 1; n < N_el - 1; ++n)
		{
			error += integrate.error(force, h1, x.segment(n*order - 1, order + 1), mesh.getTransform(n));
		}
		segment.head(order) = x.segment((N_el-1)*order-1, order);
		segment[order] = 0.0;
		error += integrate.error(force, h1, segment, mesh.getTransform(N_el-1));
	}
	else if(formulation == "Mixed")
	{
		for (int n = 0; n < N_el; ++n)
		{
			error += integrate.error(force, l2, x.segment((N_el+n)*order+1, order), mesh.getTransform(n));
		}
	}
	else if (formulation == "Mimetic")
	{
		VectorXd segment(order + 1);
		segment[0] = 0.0;
		segment.tail(order) = x.segment(N_el*order, order);
		error += integrate.error(force, h1, segment, mesh.getTransform(0));
		for (int n = 1; n < N_el - 1; ++n)
		{
			error += integrate.error(force, h1, x.segment((N_el+n)*order-1, order + 1), mesh.getTransform(n));
		}
		segment.head(order) = x.segment((2*N_el-1)*order-1, order);
		segment[order] = 0.0;
		error += integrate.error(force, h1, segment, mesh.getTransform(N_el-1));
	}

	error = sqrt(error);

	// std::cout << std::endl;
	std::cout << "Solution error: " << error << std::endl;

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

	// Also convert x to u and P vectors for the json output
	std::vector<double> u, P;
	if (formulation == "Standard")
	{
		P.resize(size);
		VectorXd::Map(&P[0], x.size()) = x;
	}
	else if(formulation == "Mixed")
	{
		P.resize(N_el*order);
		u.resize(N_el*order + 1);
		VectorXd::Map(&P[0], P.size()) = x.tail(P.size());
		VectorXd::Map(&u[0], u.size()) = x.head(u.size());
	}
	else if (formulation == "Mimetic")
	{
		P.resize(N_el*order - 1);
		u.resize(N_el*order);
		VectorXd::Map(&P[0], P.size()) = x.tail(P.size());
		VectorXd::Map(&u[0], u.size()) = x.head(u.size());
	}	
	output["u"] = u;
	output["P"] = P;
	int outargs = 4;

	if (input["Plot"]["Enable"] != "On")
	{
		std::ofstream outfile(argv[2]);
		outfile << std::setw(outargs) << output << std::endl;
	} else
	{
		// Set up from grids for plotting output
		int spe = (int)input["Plot"]["SPE"]; //samples per element
		double len = (double)input["Mesh"]["x_max"] - (double)input["Mesh"]["x_min"];
		std::vector<double> xs, xe, us, Ps, h1s, l2s;
		xs.resize(spe*N_el);
		xe.resize(spe);
		us.resize(spe*N_el);
		Ps.resize(spe*N_el);
		h1s.resize(order+1);
		l2s.resize(order);

		// Determine sample "points" per element - mesh of <spe> points from -1 to 1
		for (int i=0; i<spe; ++i)
		{
			xe[i] = 2*i*(1/(double)(spe-1))-1;
		}

		// Fill in the grids for plotting output
		Transform1D transform;
		transform = mesh.getTransform(0);
		for (int i = 0; i < spe; ++i)
		{
			xs[i] = transform.forwardTransform(xe[i]);
			VectorXd::Map(&h1s[0], h1s.size()) = h1.eval(xe[i]);
			VectorXd::Map(&l2s[0], l2s.size()) = l2.eval(xe[i]);
			if(formulation == "Mixed") us[i] += h1s[0]*u[0];
			for (int k=1; k<order+1; ++k)
			{
				if(formulation == "Mixed")
				{
					us[i] += h1s[k]*u[k];
					Ps[i] += l2s[k-1]*P[k-1];
				}
				else
				{
					Ps[i] += h1s[k]*P[k-1];
					if(formulation == "Mimetic") us[i] += l2s[k-1]*u[k-1];
				}
			}
		}
		for (int j=1; j < N_el - 1; ++j)
		{
			transform = mesh.getTransform(j);
			for (int i=0; i<spe; ++i)
			{
				xs[i+j*spe] = transform.forwardTransform(xe[i]);
				VectorXd::Map(&h1s[0], h1s.size()) = h1.eval(xe[i]);
				VectorXd::Map(&l2s[0], l2s.size()) = l2.eval(xe[i]);
				for (int k=0; k<order+1; ++k)
				{
					if(formulation == "Mixed")
					{
						us[i+j*spe] += h1s[k]*u[order*j+k];
						if(k != order) Ps[i+j*spe] += l2s[k]*P[order*j+k];
					}
					else
					{
						Ps[i+j*spe] += h1s[k]*P[order*j+(k-1)];
						if ( (formulation == "Mimetic") && (k!=order) )
						{
							us[i+j*spe] += l2s[k]*u[order*j+k];
						}
					}
				}
			}
		}
		transform = mesh.getTransform(N_el-1);
		for (int i = 0; i < spe; ++i)
		{
			xs[i+(N_el-1)*spe] = transform.forwardTransform(xe[i]);
			VectorXd::Map(&h1s[0], h1s.size()) = h1.eval(xe[i]);
			VectorXd::Map(&l2s[0], l2s.size()) = l2.eval(xe[i]);
			if(formulation == "Mixed") us[i+(N_el-1)*spe] += h1s[order]*u[order*N_el];
			for (int k=0; k<order; ++k)
			{
				if(formulation == "Mixed")
				{
					us[i+(N_el-1)*spe] += h1s[k]*u[order*(N_el-1)+k];
					Ps[i+(N_el-1)*spe] += l2s[k]*P[order*(N_el-1)+k];
				}
				else
				{
					Ps[i+(N_el-1)*spe] += h1s[k]*P[order*(N_el-1)+(k-1)];
					if(formulation == "Mimetic") us[i+(N_el-1)*spe] += l2s[k]*u[order*(N_el-1)+k];
				}
			}
		}

		output["xgrid"] = xs;
		output["Pgrid"] = Ps;
		output["ugrid"] = us;
		outargs += 3;

		std::ofstream outfile(argv[2]);
		outfile << std::setw(outargs) << output << std::endl;

		// If you want to run python plotting automatically from here
		// /*
		std::cout << std::endl;
		std::cout << "Running python plotter ..." << std::endl;
		std::string ifile = argv[1];
		std::string ofile = argv[2];
		std::string pyplot = "./plot/pyplot.py -i " + ifile + " -o " + ofile;
		std::system(pyplot.c_str());
		//*/
	}

	return 0;
}
