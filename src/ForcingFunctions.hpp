#ifndef _ForcingFunctions_hpp
#define _ForcingFunctions_hpp

class Force1D
{
	public:
		Force1D() = default;
		virtual ~Force1D() = default;
		virtual double f(double x) = 0;
};

#endif
