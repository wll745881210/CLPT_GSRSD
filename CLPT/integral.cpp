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
    x_buf.clear(  );
    y_buf.clear(  );
    return;
}

void integral::clear( const unsigned & size )
{
    x_buf.resize( size );
    y_buf.resize( size );
    return;
}

////////////////////////////////////////////////////////////
// Read data

void integral::read( const double & x, const double & y )
{
    x_buf.push_back( x );
    y_buf.push_back( y );
    return;
}

////////////////////////////////////////////////////////////
// Integrate--just in case that the x are not evenly spaced

double integral::result(  )
{
    double res( 0. );
#pragma omp parallel for reduction( + : res )
    for( unsigned i = 1; i < x_buf.size(  ); ++ i )
	res += 0.5 * ( x_buf.at( i ) - x_buf.at( i - 1 ) )
		   * ( y_buf.at( i ) + y_buf.at( i - 1 ) );
    return res;
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
    return x[ i ];
}

void integral::gl_read( const int & i, const double & kernel )
{
    gl_intg_res += w[ i ] * kernel;
    return;
}

double integral::gl_result(  )
{
    return gl_intg_res;
}
