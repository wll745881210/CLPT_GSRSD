#include "integral.h"

#include <iostream>
#include <cmath>

////////////////////////////////////////////////////////////
// Static variables

#include "gl_intg.h"

////////////////////////////////////////////////////////////
// Constructor, desctructor and initializer

integral::integral(  )
{
	clear(  );
}

integral::~integral(  )
{
	
}

void integral::clear(  )
{
	x_buf[ 0 ] = -1e-16;
	x_buf[ 1 ] = -1e-16;
	y_buf[ 0 ] = -1e-16;
	y_buf[ 1 ] = -1e-16;
	intg_res = 0.;
	counter = 0;

	return;
}

////////////////////////////////////////////////////////////
// Read data

void integral::read( const double & x, const double & y )
{
	x_buf[ 0 ] = x_buf[ 1 ];
	y_buf[ 0 ] = y_buf[ 1 ];
	x_buf[ 1 ] = x;
	y_buf[ 1 ] = y;

	const double temp_res
		= ( y_buf[ 1 ] + y_buf[ 0 ] ) * 0.5
		* fabs( x_buf[ 1 ] - x_buf[ 0 ] );
	
	intg_res += temp_res;
	
	++ counter;
	
	return;
}

////////////////////////////////////////////////////////////
// Feedback

double integral::result(  )
{
	return intg_res;
}

////////////////////////////////////////////////////////////
// Gauss-Legendre

void integral::gl_clear(  )
{
	this->gl_intg_res = 0.;
	return;
}

double integral::gl_xi( const int & i )
{
	if( i > gl_num )
		throw "Gauss-Legendre quadrature out of range";
	
	i_current = i;
	return x[ i ];
}

void integral::gl_read( const double & kernel )
{
	gl_intg_res += w[ i_current ] * kernel;	
	return;
}

double integral::gl_result(  )
{
	return gl_intg_res;
}

