#include <Eigen/Core>

#ifndef _GaussLegendre_hpp
#define _GaussLegendre_hpp

using namespace Eigen;

class GaussLegendre
{
	public:
		GaussLegendre();
		~GaussLegendre();
	private:
		MatrixXd leg2lag;
		PartialPivLU<MatrixXd> lu;
		VectorXd x, w;
};

#endif