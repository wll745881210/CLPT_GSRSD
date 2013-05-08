#ifndef INTERP1_H_
#define INTERP1_H_

#include <vector>
#include <map>

class interp1
{
	////////// Initializer //////////
public:
	interp1(  );
	~interp1(  );
	void clear(  );

	////////// Read data //////////
public:
	void read( const double & x, const double & y );
	void read( const std::vector<double> & x,
			   const std::vector<double> & y );

	////////// Interpolate //////////
private:						// Data
	std::map<double, double> xy_buf;
public:
	double operator() ( const double & arg );

	////////// Misc func/const //////////
private:
	static const double nearly0 = 1e-8;
};

#endif

