#include "Basis1D.hpp"

#ifndef _ReferenceElements_hpp
#define _ReferenceElements_hpp

class H1_1D
{
	public:
		H1_1D();
		~H1_1D();
	private:
		GaussLobatto gll;
};

#endif