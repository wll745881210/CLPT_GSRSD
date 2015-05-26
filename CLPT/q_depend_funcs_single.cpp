#include "q_depend_funcs_single.h"
#include "integral.h"
#include <iostream>
#include <cmath>
#include <map>

////////////////////////////////////////////////////////////
// Static variables

std::vector<double>   q_func_single::qvec;
std::vector<double>   q_func_single::kvec;
std::map<double, int> q_func_single::q_idx_buf;

int    q_func_single::idx_current( 0. );
double q_func_single::q_current  ( 0. );

////////////////////////////////////////////////////////////
// Constructor/destructor

q_func_single::q_func_single(  )
{
    p_kf = k_func::get_instance(  );
    return;
}

q_func_single::~q_func_single(  )
{
    return;
}

////////////////////////////////////////////////////////////
// Integration kernel is purely virtual:
// virtual double q_func_single::kernel
// ( const double & k, const double & jx,
//   const k_func & kf ) = 0;

////////////////////////////////////////////////////////////
// Integrator

double q_func_single::sph_bessel_j
( const int & n, const double & x )
{
    static const double nearly_0 = 1.e-5;
    if( fabs( x ) < nearly_0 )
	switch( n )
	{
	case 0:
	    return 1. - pow( x, 2 ) / 6.;
	case 1:
	    return x / 3. - pow( x, 3 ) / 30.;
	case 2:
	    return pow( x, 2 ) / 15.;
	case 3:
	    return pow( x, 3 ) / 105.;
	default:
	    throw "Bessel function error.";
	}
    switch( n )
    {
    case 0:
	return sin( x ) / x;
    case 1:
	return sin( x ) / pow( x, 2 ) - cos( x ) / x;
    case 2:
	return ( 3. / pow( x, 2 ) - 1 ) * sin( x ) / x
	       - 3 * cos( x ) / pow( x, 2 );
    case 3:
	return ( 15. / pow( x, 3 ) - 6./ x ) * sin( x ) / x
	     - ( 15. / pow( x, 2 ) - 1     ) * cos( x ) / x;
    default:
	throw "Bessel function error.";
    }
    return 0.;
}

double q_func_single::integrate_single( const double & q )
{
    integral intg;
    for( unsigned i = 0; i < kvec.size(  ); ++ i )
    {
	const double & k = kvec[ i ];
	const double  jx = k * q; // Argument of j_n( x )
	intg.read( k, kernel( k, jx, * p_kf ) );
    }

    static const double one_over_pi2 = 0.050660592;
    return one_over_pi2 * intg.result(  );
}

void q_func_single::integrate(  )
{
    valvec.clear(  );
    for( unsigned i = 0; i < qvec.size(  ); ++ i )
	valvec.push_back( integrate_single( qvec[ i ] ) );

    return;
}

////////////////////////////////////////////////////////////
// Integration results registery

void q_func_single::set_kvec
( const std::vector<double> & k_src )
{
    kvec = k_src;
    return;
}

void q_func_single::set_qvec
( const std::vector<double> & q_src )
{
    qvec.clear(  );
    q_idx_buf.clear(  );
    for( unsigned i = 0; i < q_src.size(  ); ++ i )
    {
	const double & q = q_src[ i ];
	const std::pair<double, int> idx_pair( q, i );
	q_idx_buf.insert( idx_pair );
	qvec.push_back( q );
    }
    return;
}

const std::vector<double> & q_func_single::get_qvec(  )
{
    return qvec;
}

double & q_func_single::valvec_at( const int & i )
{
    return valvec.at( i );
}

////////////////////////////////////////////////////////////
// Interpolation

void q_func_single::eval_all_idx( const double & q )
{
    std::map<double, int>::iterator p
	= q_idx_buf.lower_bound( q );
    if( p != q_idx_buf.begin(  ) )
	-- p;
    idx_current = p->second;
    q_current   = q;
    return;
}

double q_func_single::get_val(  )
{
    static const double nearly_0 = 1.e-3;
    
    const int    & i = idx_current;
    const double & q = q_current;
    
    if( i < 0 || i > int( qvec.size(  ) - 2 ) )
	return 0.;
    else if( ( qvec[ i + 1 ] - qvec[ i ] ) / qvec[ i ]
	       < nearly_0 )
	return valvec.at( i );
    else
	return valvec[ i ]
	    + ( valvec[ i + 1 ] - valvec[ i ] )
	    / ( qvec  [ i + 1 ] - qvec  [ i ] )
	    * ( q - qvec[ i ] );
}

