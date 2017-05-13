#include "json.hpp"
#include "Integrators.hpp"

#ifndef _Mesh2D_hpp
#define _Mesh2D_hpp

using namespace nlohmann;

class Mesh2D
{
	public:
		Mesh2D(json params);
		~Mesh2D() = default;
		Transform2D getTransform(int k);
	private:
		const int N_x_el;
		const int N_y_el;
		const int N_el;
		const int p;
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
	GaussLobatto gl((int)params["map_order"] + 1);
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
			l2_temp[j_loc*p+i_loc][0] = j_loc*p+i_loc;
		}
	}
	return l2_temp;
}

Transform2D Mesh2D::getTransform(int k)
{
	return Transform2D(transform_nodes[k]);
}

#endif // _Mesh2D_hpp
