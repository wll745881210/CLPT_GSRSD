#include "interp1.h"

#include <cmath>

////////////////////////////////////////////////////////////
// Initializer

interp1::interp1(  )
{
	clear(  );
}

interp1::~interp1(  )
{
	
}

void interp1::clear(  )
{
	xy_buf.clear(  );
	return;
}

////////////////////////////////////////////////////////////
// Read data

void interp1::read( const double & x, const double & y )
{
	std::pair<double, double> p;
	p.first  = x;
	p.second = y;
	xy_buf.insert( p );
	return;
}

void interp1::read( const std::vector<double> & x,
					const std::vector<double> & y )
{
	if( x.size(  ) != y.size(  ) )
		throw "Unable to create spline: "
			"Dimension must agree.";
	clear(  );
	for( unsigned i = 0; i < x.size(  ); ++ i )
		read( x[ i ], y[ i ] );

	return;
}

////////////////////////////////////////////////////////////
// Do interpolate

double interp1::operator() ( const double & arg )
{
	std::map<double, double>::iterator p
		= xy_buf.lower_bound( arg );
	if( p == xy_buf.end(  ) )
	{
		-- p;
		return p->second;
	}
	else if( p == xy_buf.begin(  ) )
		return p->second;
    else
	{
		const double x_1 = p->first;
		const double y_1 = p->second;
		-- p;
		const double x_0 = p->first;
		const double y_0 = p->second;
		if( fabs( x_0 - x_1 ) < nearly0 )
			return p->second;
		return y_0 + ( y_1 - y_0 )
			* ( arg - x_0 ) / ( x_1 - x_0 );
	}
}


