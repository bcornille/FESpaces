#include "json.hpp"
#include "Integrators.hpp"
#include "Eigen/Sparse"

#ifndef _Mesh2D_hpp
#define _Mesh2D_hpp

using namespace nlohmann;
using namespace Eigen;

class Mesh2D
{
	public:
		Mesh2D(json params);
		~Mesh2D() = default;
		Transform2D getTransform(int k);
		SparseMatrix<double> assembleStandard(Integrator2D integrator);
		VectorXd rhsStandard(Integrator2D integrator, const std::shared_ptr<Force2D>& f);
		double errorStandard(Integrator2D integrator, const std::shared_ptr<Force2D>& f,
			VectorXd solution);
		SparseMatrix<double> assembleMixed(Integrator2D integrator);
		VectorXd rhsMixed(Integrator2D integrator, const std::shared_ptr<Force2D>& f);
		double errorMixed(Integrator2D integrator, const std::shared_ptr<Force2D>& f,
			VectorXd solution);
		SparseMatrix<double> assembleDualMixed(Integrator2D integrator);
		VectorXd rhsDualMixed(Integrator2D integrator, const std::shared_ptr<Force2D>& f);
		double errorDualMixed(Integrator2D integrator, const std::shared_ptr<Force2D>& f,
			VectorXd solution);
		auto  samplePStandard(VectorXd solution, json params, json plot);
		auto  samplePuDualMixed(VectorXd solution, json params, json plot);
		auto  samplePuMixed(VectorXd solution, json params, json plot);
	private:
		const int N_x_el;
		const int N_y_el;
		const int N_el;
		const int p;
		const int p_map;
		std::vector<Matrix2Xd> transform_nodes;
		std::vector<std::vector<Vector2i> > h1_dofs;
		std::vector<std::vector<Vector2i> > hcurl_dofs;
		std::vector<std::vector<Vector2i> > hdiv_dofs;
		std::vector<std::vector<Vector2i> > l2_dofs;
		H1_2D h1_el;
		HCurl_2D hcurl_el;
		HDiv_2D hdiv_el;
		L2_2D l2_el;
		int N_h1_dofs;
		int N_hcurl_dofs;
		int N_hdiv_dofs;
		int N_l2_dofs;
		std::vector<Vector2i> h1_dofs_setup(int i, int j);
		std::vector<Vector2i> hcurl_dofs_setup(int i, int j);
		std::vector<Vector2i> hdiv_dofs_setup(int i, int j);
		std::vector<Vector2i> l2_dofs_setup(int i, int j);
};

Mesh2D::Mesh2D(json params) :
	N_x_el((int)params["N_x"]), N_y_el((int)params["N_y"]),
	N_el((int)params["N_x"]*(int)params["N_y"]), p((int)params["el_order"]),
	p_map((int)params["map_order"]),
	transform_nodes((int)params["N_x"]*(int)params["N_y"]),
	h1_dofs((int)params["N_x"]*(int)params["N_y"]),
	hcurl_dofs((int)params["N_x"]*(int)params["N_y"]),
	hdiv_dofs((int)params["N_x"]*(int)params["N_y"]),
	l2_dofs((int)params["N_x"]*(int)params["N_y"]),
	h1_el((int)params["el_order"]), hcurl_el((int)params["el_order"]),
	hdiv_el((int)params["el_order"]), l2_el((int)params["el_order"])
{
	double x_len = (double)params["x_max"] - (double)params["x_min"];
	double y_len = (double)params["y_max"] - (double)params["y_min"];
	double delta_xi = 2.0/N_x_el;
	double delta_eta = 2.0/N_y_el;
	double c = (double)params["c"];
	GaussLobatto gl(p_map + 1);
	int N_map = gl.getN();
	Matrix2Xd temp_nodes(2, N_map*N_map);
	N_h1_dofs = 0;
	N_hcurl_dofs = 0;
	N_hdiv_dofs = 0;
	N_l2_dofs = 0;
	for (int j = 0; j < N_y_el; ++j)
	{
		for (int i = 0; i < N_x_el; ++i)
		{
			int k = j*N_x_el + i;
			h1_dofs[k] = h1_dofs_setup(i, j);
			hcurl_dofs[k] = hcurl_dofs_setup(i, j);
			hdiv_dofs[k] = hdiv_dofs_setup(i, j);
			l2_dofs[k] = l2_dofs_setup(i, j);
			Transform1D x_segment(-1.0 + i*delta_xi, -1.0 + (i + 1)*delta_xi);
			Transform1D y_segment(-1.0 + j*delta_eta, -1.0 + (j + 1)*delta_eta);
			for (int j_loc = 0; j_loc < N_map; ++j_loc)
			{
				double eta = y_segment.forwardTransform(gl.getNode(j_loc));
				for (int i_loc = 0; i_loc < N_map; ++i_loc)
				{
					double xi = x_segment.forwardTransform(gl.getNode(i_loc));
					temp_nodes(0, j_loc*N_map + i_loc) = x_len*(1.0 + xi
						+ c*sin(pi()*xi)*sin(pi()*eta))/2.0;
					temp_nodes(1, j_loc*N_map + i_loc) = y_len*(1.0 + eta
						+ c*sin(pi()*xi)*sin(pi()*eta))/2.0;
				}
			}
			transform_nodes[k] = temp_nodes;
		}
	}
}

std::vector<Vector2i> Mesh2D::h1_dofs_setup(int i, int j)
{
	std::vector<Vector2i> h1_temp;
	int k = j*N_x_el + i;
	int i_loc_start = 0;
	if (i == 0) i_loc_start = 1;
	int j_loc_start = 0;
	if (j == 0) j_loc_start = 1;
	int i_loc_end = p;
	if (i == N_x_el - 1) i_loc_end = p - 1;
	int j_loc_end = p;
	if (j == N_y_el - 1) j_loc_end = p - 1;
	int n_i = i_loc_end - i_loc_start + 1;
	h1_temp.resize(n_i*(j_loc_end - j_loc_start + 1));
	for (int j_loc = j_loc_start; j_loc < j_loc_end + 1; ++j_loc)
	{
		int j_offset = j_loc - j_loc_start;
		for (int i_loc = i_loc_start; i_loc < i_loc_end + 1; ++i_loc)
		{
			if (j_loc == 0)
			{
				if (i_loc == 0)
				{
					if (i == 1)
					{
						h1_temp[0][0] = h1_dofs[k-1][p-1][0];
					}
					else
					{
						h1_temp[0][0] = h1_dofs[k-1][p][0];
					}
				}
				else
				{
					if (i == 0)
					{
						if (j == 1)
						{
							if (i == N_x_el - 1)
							{
								h1_temp[i_loc-i_loc_start][0] = h1_dofs[k-N_x_el][(p-1)*(p-1)+i_loc-1][0];
							}
							else
							{
								h1_temp[i_loc-i_loc_start][0] = h1_dofs[k-N_x_el][p*(p-1)+i_loc-1][0];
							}
						}
						else
						{
							if (i == N_x_el - 1)
							{
								h1_temp[i_loc-i_loc_start][0] = h1_dofs[k-N_x_el][(p-1)*p+i_loc-1][0];
							}
							else
							{
								h1_temp[i_loc-i_loc_start][0] = h1_dofs[k-N_x_el][p*p+i_loc-1][0];
							}
						}
					}
					else if (j == 1)
					{
						if (i == N_x_el - 1)
						{
							h1_temp[i_loc-i_loc_start][0] = h1_dofs[k-N_x_el][p*(p-1)+i_loc][0];
						}
						else
						{
							h1_temp[i_loc-i_loc_start][0] = h1_dofs[k-N_x_el][(p+1)*(p-1)+i_loc][0];
						}
					}
					else
					{
						if (i == N_x_el - 1)
						{
							h1_temp[i_loc-i_loc_start][0] = h1_dofs[k-N_x_el][p*p+i_loc][0];
						}
						else
						{
							h1_temp[i_loc-i_loc_start][0] = h1_dofs[k-N_x_el][(p+1)*p+i_loc][0];
						}
					}
				}
			}
			else if (i_loc == 0)
			{
				if (j == 0)
				{
					if (i == 1)
					{
						h1_temp[j_offset*n_i][0] = h1_dofs[k-1][p*j_loc-1][0];
					}
					else
					{
						h1_temp[j_offset*n_i][0] = h1_dofs[k-1][(p+1)*j_loc-1][0];
					}
				}
				else if (i == 1)
				{
					h1_temp[j_offset*n_i][0] = h1_dofs[k-1][p*(j_loc+1)-1][0];
				}
				else
				{
					h1_temp[j_offset*n_i][0] = h1_dofs[k-1][(p+1)*(j_loc+1)-1][0];
				}
			}
			else
			{
				h1_temp[j_offset*n_i+i_loc-i_loc_start][0] = N_h1_dofs++;
			}
			h1_temp[j_offset*n_i+i_loc-i_loc_start][1] = j_loc*(p + 1) + i_loc;
		}
	}
	return h1_temp;
}

std::vector<Vector2i> Mesh2D::hcurl_dofs_setup(int i, int j)
{
	std::vector<Vector2i> hcurl_temp;
	int k = j*N_x_el + i;
	hcurl_temp.resize(2*p*(p + 1));
	for (int j_loc = 0; j_loc < p + 1; ++j_loc)
	{
		for (int i_loc = 0; i_loc < p; ++i_loc)
		{
			if (j_loc == 0)
			{
				if (j == 0)
				{
					hcurl_temp[j_loc*p+i_loc][0] = N_hcurl_dofs++;
				}
				else
				{
					hcurl_temp[j_loc*p+i_loc][0] = hcurl_dofs[k-N_x_el][p*p+i_loc][0];
				}
			}
			else
			{
				hcurl_temp[j_loc*p+i_loc][0] = N_hcurl_dofs++;
			}
			hcurl_temp[j_loc*p+i_loc][1] = j_loc*p + i_loc;
		}
	}
	for (int j_loc = 0; j_loc < p; ++j_loc)
	{
		for (int i_loc = 0; i_loc < p + 1; ++i_loc)
		{
			if (i_loc == 0)
			{
				if (i == 0)
				{
					hcurl_temp[(j_loc+p)*(p+1)+i_loc][0] = N_hcurl_dofs++;
				}
				else
				{
					hcurl_temp[(j_loc+p)*(p+1)+i_loc][0] = hcurl_dofs[k-1][(p+1)*(p+j_loc+1)-1][0];
				}
			}
			else
			{
				hcurl_temp[(j_loc+p)*(p+1)+i_loc][0] = N_hcurl_dofs++;
			}
			hcurl_temp[(j_loc+p)*(p+1)+i_loc][1] = (p + j_loc)*(p + 1) + i_loc;
		}
	}
	return hcurl_temp;
}

std::vector<Vector2i> Mesh2D::hdiv_dofs_setup(int i, int j)
{
	std::vector<Vector2i> hdiv_temp;
	int k = j*N_x_el + i;
	hdiv_temp.resize(2*p*(p + 1));
	for (int j_loc = 0; j_loc < p + 1; ++j_loc)
	{
		for (int i_loc = 0; i_loc < p; ++i_loc)
		{
			if (j_loc == 0)
			{
				if (j == 0)
				{
					hdiv_temp[j_loc*p+i_loc][0] = N_hdiv_dofs++;
				}
				else
				{
					hdiv_temp[j_loc*p+i_loc][0] = hdiv_dofs[k-N_x_el][p*p+i_loc][0];
				}
			}
			else
			{
				hdiv_temp[j_loc*p+i_loc][0] = N_hdiv_dofs++;
			}
			hdiv_temp[j_loc*p+i_loc][1] = j_loc*p + i_loc;
		}
	}
	for (int j_loc = 0; j_loc < p; ++j_loc)
	{
		for (int i_loc = 0; i_loc < p + 1; ++i_loc)
		{
			if (i_loc == 0)
			{
				if (i == 0)
				{
					hdiv_temp[(j_loc+p)*(p+1)+i_loc][0] = N_hdiv_dofs++;
				}
				else
				{
					hdiv_temp[(j_loc+p)*(p+1)+i_loc][0] = hdiv_dofs[k-1][(p+1)*(p+j_loc+1)-1][0];
				}
			}
			else
			{
				hdiv_temp[(j_loc+p)*(p+1)+i_loc][0] = N_hdiv_dofs++;
			}
			hdiv_temp[(j_loc+p)*(p+1)+i_loc][1] = (p + j_loc)*(p + 1) + i_loc;
		}
	}
	return hdiv_temp;
}

std::vector<Vector2i> Mesh2D::l2_dofs_setup(int i, int j)
{
	std::vector<Vector2i> l2_temp;
	l2_temp.resize(p*p);
	for (int j_loc = 0; j_loc < p; ++j_loc)
	{
		for (int i_loc = 0; i_loc < p; ++i_loc)
		{
			l2_temp[j_loc*p+i_loc][0] = N_l2_dofs++;
			l2_temp[j_loc*p+i_loc][1] = j_loc*p+i_loc;
		}
	}
	return l2_temp;
}

Transform2D Mesh2D::getTransform(int k)
{
	return Transform2D(transform_nodes[k], p_map);
}

SparseMatrix<double> Mesh2D::assembleStandard(Integrator2D integrator)
{
	SparseMatrix<double> matrix(N_h1_dofs, N_h1_dofs);
	matrix.reserve(VectorXi::Constant(N_h1_dofs, (2*p + 1)*(2*p + 1)));
	matrix.setZero();
	for (int k = 0; k < N_el; ++k)
	{
		MatrixXd minimatrix = integrator.laplace(h1_el, h1_el, getTransform(k));
		for (Vector2i& dof_j : h1_dofs[k])
		{
			for (Vector2i& dof_i : h1_dofs[k])
			{
				matrix.coeffRef(dof_i[0], dof_j[0]) += minimatrix(dof_i[1], dof_j[1]);
			}
		}
	}
	matrix.makeCompressed();
	return matrix;
}

VectorXd Mesh2D::rhsStandard(Integrator2D integrator, const std::shared_ptr<Force2D>& f)
{
	VectorXd rhs(N_h1_dofs);
	rhs.setZero();
	for (int k = 0; k < N_el; ++k)
	{
		VectorXd minirhs = integrator.force(f, h1_el, getTransform(k));
		for (Vector2i& dof_i : h1_dofs[k])
		{
			rhs[dof_i[0]] += minirhs[dof_i[1]];
		}
	}
	return rhs;
}

double Mesh2D::errorStandard(Integrator2D integrator, const std::shared_ptr<Force2D>& f,
	VectorXd solution)
{
	double error = 0.0;
	VectorXd coeffs((p + 1)*(p + 1));
	for (int k = 0; k < N_el; ++k)
	{
		coeffs.setZero();
		for (Vector2i& dof_i : h1_dofs[k])
		{
			coeffs[dof_i[1]] = solution[dof_i[0]];
		}
		error += integrator.error(f, h1_el, coeffs, getTransform(k));
	}
	return error;
}

auto Mesh2D::samplePStandard(VectorXd solution, json params, json plot)
{
	// Return a surface defined by 2D arrays x, y, and P
	struct surface 
	{
		std::vector<std::vector<double> > x; 
		std::vector<std::vector<double> > y; 
		std::vector<std::vector<double> > P;
	};

	// Set up the grid
	const int spef = (int)plot["SPEf"];
	const int Npef = spef*spef; 
	std::vector<std::vector<double> > xgrid(spef*N_x_el, std::vector<double>(spef*N_y_el));
	std::vector<std::vector<double> > ygrid(spef*N_x_el, std::vector<double>(spef*N_y_el));
	std::vector<std::vector<double> > Pgrid(spef*N_x_el, std::vector<double>(spef*N_y_el));
	std::vector<std::vector<Vector2d> > xys(spef, std::vector<Vector2d>(spef));
	double dX = ((double)params["x_max"] - (double)params["x_min"])/(N_x_el);
	double dY = ((double)params["y_max"] - (double)params["y_min"])/(N_y_el);
	double dx = dX/(spef-1);
	double dy = dY/(spef-1);
	VectorXd coeffs(h1_el.dofs());
	Vector2d loc;
	Transform2D transform; 

	// Determine sample "points" per element - mesh of <spe> points from -1 to 1
	for (int m = 0; m < spef; ++m)
	{
		for (int n = 0; n < spef; ++n)
		{
			xys[n][m][0] = 2*n*(1.0/(double)(spef-1.0)) - 1;
			xys[n][m][1] = 2*m*(1.0/(double)(spef-1.0)) - 1;
		}
	}

	// Loop through cells
	for (int i = 0; i < N_y_el; ++i)
	{
		for (int j = 0; j < N_x_el; ++j)
		{
			// Determine coefficients for the cell we are in
			int k = j + N_x_el*i;
			coeffs.setZero();
			for (Vector2i& dof_i : h1_dofs[k])
			{
				coeffs[dof_i[1]] = solution[dof_i[0]];
			}
			
			transform = getTransform(k);

			// Determine x, y, and P values for the cell we are in
			for (int m = 0; m < spef; ++m)
			{
				for (int n = 0; n < spef; ++n)
				{
					Vector2d xystemp = transform.forwardTransform(xys[n][m]);
					xgrid[n+j*spef][m+i*spef] = xystemp[0];
					ygrid[n+j*spef][m+i*spef] = xystemp[1];
					loc = xys[n][m];
					Pgrid[n+j*spef][m+i*spef] = coeffs.dot(h1_el.eval(loc));
				}
			}
		}
	}

	return surface {xgrid, ygrid, Pgrid};
}

SparseMatrix<double> Mesh2D::assembleMixed(Integrator2D integrator)
{
	SparseMatrix<double> matrix(N_hdiv_dofs + N_l2_dofs, N_hdiv_dofs + N_l2_dofs);
	matrix.reserve(VectorXi::Constant(N_hdiv_dofs + N_l2_dofs, 3*p*(2*p + 1)));
	matrix.setZero();
	for (int k = 0; k < N_el; ++k)
	{
		MatrixXd mini_hdiv_mass = integrator.mass(hdiv_el, hdiv_el, getTransform(k));
		MatrixXd mini_div = integrator.div(hdiv_el, l2_el, getTransform(k));
		for (Vector2i& dof_j : hdiv_dofs[k])
		{
			for (Vector2i& dof_i : hdiv_dofs[k])
			{
				matrix.coeffRef(dof_i[0], dof_j[0]) -= mini_hdiv_mass(dof_i[1], dof_j[1]);
			}
			for (Vector2i& dof_i : l2_dofs[k])
			{
				matrix.coeffRef(N_hdiv_dofs + dof_i[0], dof_j[0]) += mini_div(dof_i[1], dof_j[1]);
			}
		}
	}
	SparseMatrix<double> identity(N_hdiv_dofs + N_l2_dofs, N_hdiv_dofs + N_l2_dofs);
	identity.setIdentity();
	matrix = matrix.selfadjointView<Lower>()*identity;
	matrix.makeCompressed();
	return matrix;
}

VectorXd Mesh2D::rhsMixed(Integrator2D integrator, const std::shared_ptr<Force2D>& f)
{
	VectorXd rhs(N_hdiv_dofs + N_l2_dofs);
	rhs.setZero();
	for (int k = 0; k < N_el; ++k)
	{
		VectorXd minirhs = integrator.force(f, l2_el, getTransform(k));
		for (Vector2i& dof_i : l2_dofs[k])
		{
			rhs[N_hdiv_dofs+dof_i[0]] += minirhs[dof_i[1]];
		}
	}
	return rhs;
}

double Mesh2D::errorMixed(Integrator2D integrator, const std::shared_ptr<Force2D>& f,
	VectorXd solution)
{
	double error = 0.0;
	VectorXd coeffs(p*p);
	for (int k = 0; k < N_el; ++k)
	{
		coeffs.setZero();
		for (Vector2i& dof_i : l2_dofs[k])
		{
			coeffs[dof_i[1]] = solution[N_hdiv_dofs+dof_i[0]];
		}
		error += integrator.error(f, l2_el, coeffs, getTransform(k));
	}
	return error;
}

auto Mesh2D::samplePuMixed(VectorXd solution, json params, json plot)
{
	// Return a surface defined by 2D arrays x, y, ux, uy, P
	struct surface 
	{
		std::vector<std::vector<double> > x; 
		std::vector<std::vector<double> > y; 
		std::vector<std::vector<double> > P;
		std::vector<std::vector<double> > ux;
		std::vector<std::vector<double> > uy;
	};

	// Set up the grid
	const int spef = (int)plot["SPEf"];
	const int Npef = spef*spef; 
	std::vector<std::vector<double> > xgrid(spef*N_x_el, std::vector<double>(spef*N_y_el));
	std::vector<std::vector<double> > ygrid(spef*N_x_el, std::vector<double>(spef*N_y_el));
	std::vector<std::vector<double> > Pgrid(spef*N_x_el, std::vector<double>(spef*N_y_el));
	std::vector<std::vector<double> > uxgrid(spef*N_x_el, std::vector<double>(spef*N_y_el));
	std::vector<std::vector<double> > uygrid(spef*N_x_el, std::vector<double>(spef*N_y_el));
	std::vector<std::vector<Vector2d> > xys(spef, std::vector<Vector2d>(spef));
	double dX = ((double)params["x_max"] - (double)params["x_min"])/(N_x_el);
	double dY = ((double)params["y_max"] - (double)params["y_min"])/(N_y_el);
	double dx = dX/(spef-1);
	double dy = dY/(spef-1);
	VectorXd coeffsP(l2_el.dofs());
	RowVectorXd coeffsu(hdiv_el.dofs());
	Vector2d loc;
	Transform2D transform;

	// Determine sample "points" per element - mesh of <spe> points from -1 to 1
	for (int m = 0; m < spef; ++m)
	{
		for (int n = 0; n < spef; ++n)
		{
			xys[n][m][0] = 2*n*(1.0/(double)(spef-1.0)) - 1;
			xys[n][m][1] = 2*m*(1.0/(double)(spef-1.0)) - 1;
		}
	}

	// Loop through cells
	for (int i = 0; i < N_y_el; ++i)
	{
		for (int j = 0; j < N_x_el; ++j)
		{
			// Determine coefficients for the cell we are in
			int k = j + N_x_el*i;
			coeffsP.setZero();
			coeffsu.setZero();
			for (Vector2i& dof_i : l2_dofs[k])
			{
				coeffsP[dof_i[1]] = solution[N_hdiv_dofs+dof_i[0]];
			}
			for (Vector2i& dof_i : hdiv_dofs[k])
			{
				coeffsu[dof_i[1]] = solution[dof_i[0]];
			}
						
			transform = getTransform(k);

			// Determine x, y, and P values for the cell we are in
			for (int m = 0; m < spef; ++m)
			{
				for (int n = 0; n < spef; ++n)
				{
					Vector2d xystemp = transform.forwardTransform(xys[n][m]);
					xgrid[n+j*spef][m+i*spef] = xystemp[0];
					ygrid[n+j*spef][m+i*spef] = xystemp[1];
					loc = xys[n][m];
					Pgrid[n+j*spef][m+i*spef] = coeffsP.dot(l2_el.eval(loc));

					Vector2d utemp = ((coeffsu*hdiv_el.eval(loc)) * 
						(transform.jacobianMatrix(xys[n][m])).transpose()).transpose();

					uxgrid[n+j*spef][m+i*spef] = utemp[0]/(transform.jacobian(xys[n][m]));
					uygrid[n+j*spef][m+i*spef] = utemp[1]/(transform.jacobian(xys[n][m]));
				}
			}
		}
	}

	return surface {xgrid, ygrid, Pgrid, uxgrid, uygrid};
}

SparseMatrix<double> Mesh2D::assembleDualMixed(Integrator2D integrator)
{
	SparseMatrix<double> matrix(N_hcurl_dofs + N_h1_dofs, N_hcurl_dofs + N_h1_dofs);
	matrix.reserve(VectorXi::Constant(N_hcurl_dofs + N_h1_dofs, 6*p*(p + 1) + 1));
	matrix.setZero();
	for (int k = 0; k < N_el; ++k)
	{
		MatrixXd mini_hcurl_mass = integrator.mass(hcurl_el, hcurl_el, getTransform(k));
		MatrixXd mini_div = integrator.div(hcurl_el, h1_el, getTransform(k));
		for (Vector2i& dof_j : hcurl_dofs[k])
		{
			for (Vector2i& dof_i : hcurl_dofs[k])
			{
				matrix.coeffRef(dof_i[0], dof_j[0]) += mini_hcurl_mass(dof_i[1], dof_j[1]);
			}
			for (Vector2i& dof_i : h1_dofs[k])
			{
				matrix.coeffRef(N_hcurl_dofs + dof_i[0], dof_j[0]) += mini_div(dof_i[1], dof_j[1]);
			}
		}
	}
	SparseMatrix<double> identity(N_hcurl_dofs + N_h1_dofs, N_hcurl_dofs + N_h1_dofs);
	identity.setIdentity();
	matrix = matrix.selfadjointView<Lower>()*identity;
	matrix.makeCompressed();
	return matrix;
}

VectorXd Mesh2D::rhsDualMixed(Integrator2D integrator, const std::shared_ptr<Force2D>& f)
{
	VectorXd rhs(N_hcurl_dofs + N_h1_dofs);
	rhs.setZero();
	for (int k = 0; k < N_el; ++k)
	{
		VectorXd minirhs = integrator.force(f, h1_el, getTransform(k));
		for (Vector2i& dof_i : h1_dofs[k])
		{
			rhs[N_hcurl_dofs+dof_i[0]] -= minirhs[dof_i[1]];
		}
	}
	return rhs;
}

double Mesh2D::errorDualMixed(Integrator2D integrator, const std::shared_ptr<Force2D>& f,
	VectorXd solution)
{
	double error = 0.0;
	VectorXd coeffs((p + 1)*(p + 1));
	for (int k = 0; k < N_el; ++k)
	{
		coeffs.setZero();
		for (Vector2i& dof_i : h1_dofs[k])
		{
			coeffs[dof_i[1]] = solution[N_hcurl_dofs+dof_i[0]];
		}
		error += integrator.error(f, h1_el, coeffs, getTransform(k));
	}
	return error;
}

auto Mesh2D::samplePuDualMixed(VectorXd solution, json params, json plot)
{
	// Return a surface defined by 2D arrays x, y, ux, uy, P
	struct surface 
	{
		std::vector<std::vector<double> > x; 
		std::vector<std::vector<double> > y; 
		std::vector<std::vector<double> > P;
		std::vector<std::vector<double> > ux;
		std::vector<std::vector<double> > uy;
	};

	// Set up the grid
	const int spef = (int)plot["SPEf"];
	const int Npef = spef*spef; 
	std::vector<std::vector<double> > xgrid(spef*N_x_el, std::vector<double>(spef*N_y_el));
	std::vector<std::vector<double> > ygrid(spef*N_x_el, std::vector<double>(spef*N_y_el));
	std::vector<std::vector<double> > Pgrid(spef*N_x_el, std::vector<double>(spef*N_y_el));
	std::vector<std::vector<double> > uxgrid(spef*N_x_el, std::vector<double>(spef*N_y_el));
	std::vector<std::vector<double> > uygrid(spef*N_x_el, std::vector<double>(spef*N_y_el));
	std::vector<std::vector<Vector2d> > xys(spef, std::vector<Vector2d>(spef));
	double dX = ((double)params["x_max"] - (double)params["x_min"])/(N_x_el);
	double dY = ((double)params["y_max"] - (double)params["y_min"])/(N_y_el);
	double dx = dX/(spef-1);
	double dy = dY/(spef-1);
	VectorXd coeffsP(h1_el.dofs());
	RowVectorXd coeffsu(hcurl_el.dofs());
	Vector2d loc;
	Transform2D transform;

	// Determine sample "points" per element - mesh of <spe> points from -1 to 1
	for (int m = 0; m < spef; ++m)
	{
		for (int n = 0; n < spef; ++n)
		{
			xys[n][m][0] = 2*n*(1.0/(double)(spef-1.0)) - 1;
			xys[n][m][1] = 2*m*(1.0/(double)(spef-1.0)) - 1;
		}
	}

	// Loop through cells
	for (int i = 0; i < N_y_el; ++i)
	{
		for (int j = 0; j < N_x_el; ++j)
		{
			// Determine coefficients for the cell we are in
			int k = j + N_x_el*i;
			coeffsP.setZero();
			coeffsu.setZero();
			for (Vector2i& dof_i : h1_dofs[k])
			{
				coeffsP[dof_i[1]] = solution[N_hcurl_dofs+dof_i[0]];
			}
			for (Vector2i& dof_i : hcurl_dofs[k])
			{
				coeffsu[dof_i[1]] = solution[dof_i[0]];
			}
						
			transform = getTransform(k);

			// Determine x, y, and P values for the cell we are in
			for (int m = 0; m < spef; ++m)
			{
				for (int n = 0; n < spef; ++n)
				{
					Vector2d xystemp = transform.forwardTransform(xys[n][m]);
					xgrid[n+j*spef][m+i*spef] = xystemp[0];
					ygrid[n+j*spef][m+i*spef] = xystemp[1];
					loc = xys[n][m];
					Pgrid[n+j*spef][m+i*spef] = coeffsP.dot(h1_el.eval(loc));

					Vector2d utemp = ((coeffsu*hcurl_el.eval(loc)) * 
						transform.jacobianInverse(xys[n][m])).transpose();

					uxgrid[n+j*spef][m+i*spef] = utemp[0];
					uygrid[n+j*spef][m+i*spef] = utemp[1];
				}
			}
		}
	}

	return surface {xgrid, ygrid, Pgrid, uxgrid, uygrid};
}

#endif // _Mesh2D_hpp
