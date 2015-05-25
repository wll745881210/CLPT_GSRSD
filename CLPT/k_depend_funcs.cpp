#include "k_depend_funcs.h"

#include <iostream>
#include <fstream>
#include <cmath>

////////////////////////////////////////////////////////////
// Static variables

const double	k_func::nearly_0     = 1.e-3;
const double	k_func::nearly_inf   = 1.e2;
const double	k_func::pi	     = 3.14159265359;
const double	k_func::one_over_pi2 = 0.050660592;

////////////////////////////////////////////////////////////
// Constructor and desctructor

k_func::k_func(  )
{

}

k_func::~k_func(  )
{
	
}

////////////////////////////////////////////////////////////
// Interpolation

double k_func::interp( const double & k,
		       const std::vector<double> & vec )
{
    itr_idx_map p = idx_map.lower_bound( k );
    if( p != idx_map.begin(  ) )
	-- p;
	
    const unsigned i = p->second;
    if( i > k_buf.size(  ) - 3 || i < 1 )
	return 0.;
	
    if( vec[ i ] > 0 && vec[ i+1 ] > 0)
    {
	const double logres
	    = log( vec[ i ] )
	    + ( log( vec[ i+1 ] ) - log( vec[ i ] ) )
	    / ( log( k_buf[ i+1 ] ) - log( k_buf[ i ] ) )
	    * ( log( k ) - log( k_buf[ i ] ) );
	return exp( logres );
    }
    else
	return vec[ i ]	+ ( vec[ i+1 ] - vec[ i ] )
	    / ( k_buf[ i+1 ] - k_buf[ i ] )
	    * ( k - k_buf[ i ] );
}

////////////////////////////////////////////////////////////
// Linear power spectrum

const std::vector<double> & k_func::kvec(  )
{
    return k_buf;
}

void k_func::load_PL( std::string file_name )
{
    std::ifstream fin( file_name.c_str(  ) );
    if( !fin )
	throw "Unable to open PL file.";

    double k( 1e-5 ), PL( 1e-5 );	
    std::pair<double, int> k_idx;
    k_idx.first  = 1e-5;
    k_idx.second = 0;

    while( !fin.eof(  ) )
    {
	k_buf.push_back( k );
	PL_buf.push_back( PL );
	idx_map.insert( k_idx );
	
	fin >> k >> PL;
	k_idx.first  = k;
	k_idx.second = k_buf.size(  );
    }
	
    k_min = k_buf[ 0 ];
    k_max = k_buf[ k_buf.size(  ) - 1 ];
    
    std::cout << "Linear power spectrum loaded from: "
	      << file_name << std::endl;
    return;
}

double k_func::PL_val( const double & k )
{
    if( k < k_min || k > k_max )
	return 0.;
    return interp( k, PL_buf );
}

////////////////////////////////////////////////////////////
// Q functions

double k_func::Q_val( const int n, const double & k )
{
    switch( n )
    {
    case 1:
	return interp( k, Q_buf_1 );
    case 2:
	return interp( k, Q_buf_2 );
    case 3:
	return interp( k, Q_buf_3 );
    case 4:
	return interp( k, Q_buf_4 );
    case 5:
	return interp( k, Q_buf_5 );
    case 6:
	return interp( k, Q_buf_6 );
    case 7:
	return interp( k, Q_buf_7 );
    case 8:
	return interp( k, Q_buf_8 );
    default:
	throw "Unable to evaluate Q function.";
    }
    return 0.;
}

void k_func::get_Q_func(  )
{
    dvec * p[ 9 ] 
	= { &Q_buf_1, &Q_buf_1, &Q_buf_2, 
	    &Q_buf_3, &Q_buf_4, &Q_buf_5,
	    &Q_buf_6, &Q_buf_7, &Q_buf_8 };
    for( int i = 0; i < 9; ++ i )
	p[ i ]->clear(  );

    std::cout << "Generating Q_n ( 1 to 8 )... "
	      << std::flush;
    // pg.init( k_buf.size(  ) );
#pragma omp parallel for
    for( int n = 1; n < 9; ++ n )
    {
	// pg.show( i );
	std::vector<double> & pn = *( p[ n ] );
	for( unsigned i = 0; i < k_buf.size(  ); ++ i )
	{
	    const double & k = k_buf[ i ];
	    const double  Qn = Q_outer_integration( n, k );
	    pn.push_back( Qn );
	}		
    }

    std::cout << "Done." << std::endl;
    return;
}

double k_func::Q_kernel
( int n, const double & r, const double & x )
{	
    if( fabs( r - 1 ) < nearly_0 )
	switch( n )
	{
	case 1:
	    return 0.25 * pow( x + 1., 2 );
	case 2:
	    return 0.25 * x * ( x + 1. );
	case 3:
	    return 0.25 * pow( x, 2 );
	case 4:
	    return ( fabs( x - 1. ) < nearly_0 ? 0. 
		     : - 0.25 * ( x + 1. ) / ( x - 1. ) );
	case 5:
	    return 0.5 * ( x + 1. ) * x	;
	case 6:
	    return 0.5 - x - 1.5 * pow( x, 2 );
	case 7:
	    return 0.5 * pow( x, 2 );
	case 8:
	    return 0.5 * ( x + 1. );
	default:
	    throw "Q kernel error.";
	}
    else if( fabs( r + 1 ) < nearly_0 )
	switch( n )
	{
	case 1:
	    return 0.25 * pow( x - 1., 2 );
	case 2:
	    return 0.25 * x * ( x - 1. );
	case 3:
	    return 0.25 * pow( x, 2 );
	case 4:
	    return ( fabs( x + 1 ) < nearly_0 ? 0. 
		     : - 0.25 * ( x - 1. ) / ( x + 1. ) );
	case 5:
	    return 0.5 * ( x - 1. ) * x;
	case 6:
	    return 0.5 + x - 1.5 * pow( x, 2 );
	case 7:
	    return 0.5 * pow( x, 2 );
	case 8:
	    return -0.5 * ( x - 1. );		
	default:
	    throw "Q kernel error.";	
	}
    else
    {
	const double y = 1. + pow( r, 2 ) - 2. * r * x;
	switch( n )
	{
	case 1:
	    return pow( r, 2 ) * pow( 1 - x*x, 2 )
		/ pow( y, 2 );
	case 2:
	    return ( 1. - x*x ) * r * x * ( 1 - r * x )
		/ pow( y, 2 );
	case 3:
	    return pow( x, 2 ) * pow( 1 - r * x, 2 )
		/ pow( y, 2 );
	case 4:
	    return ( 1. - x*x ) / pow( y, 2 );
	case 5:
	    return r * x * ( 1. - x*x ) / y;
	case 6:
	    return ( 1. - 3. * r * x ) * ( 1. - x*x ) / y;
	case 7:
	    return pow( x, 2 ) * ( 1 - r * x ) / y;
	case 8:
	    return pow( r, 2 ) * ( 1 - x*x ) / y;	
	default:
	    throw "Q kernel error.";
	}
    }
    return 0.;
}

double k_func::Q_inner_integration
( int n, const double & r, const double & k )
{
    integral intg;
    intg.gl_clear(  );
    for( int i = 0; i < intg.gl_num; ++ i )
    {
	const double x_i = intg.gl_xi( i );
	const double y   = sqrt( 1. + pow( r, 2 )
	                         - 2. * r * x_i ) * k;
	intg.gl_read( i, Q_kernel( n, r, x_i )
	                 * PL_val( y ) );
    }

    return intg.gl_result(  );	
}

double k_func::Q_outer_integration
( int n, const double & k )
{
    integral intg;
    intg.clear(  );
    for( unsigned i = 0; i < k_buf.size(  ); ++ i )
    {
	const double r      = k_buf[ i ] / k;
	const double kernel = PL_buf[ i ]
	    * Q_inner_integration( n, r, k );
	intg.read( r, kernel );
    }
	
    static const double temp_coef = 0.025330295910584444;
    // 0.25 / \pi^2;
    return pow( k, 3 ) * temp_coef * intg.result(  );
}

////////////////////////////////////////////////////////////
// R functions

double k_func::R_val( const int n, const double & k )
{
    switch( n )
    {
    case 1:
	return interp( k, R_buf_1 );
    case 2:
	return interp( k, R_buf_2 );
    default:
	throw "Unable to evaluate Q function.";
    }
    return 0.;
}

void k_func::get_R_func(  )
{
    dvec * p[ 3 ] = { &R_buf_1, &R_buf_1, &R_buf_2 };
    for( int i = 0; i < 3; ++ i )
	p[ i ]->clear(  );

    std::cout << "Generating R_n ( 1 to 2 )... "
	      << std::flush;
    // pg.init( k_buf.size(  ) );

#pragma omp parallel for
    for( int n = 1; n <= 2; ++ n )
    {
	// pg.show( i );
	for( unsigned i = 0; i < k_buf.size(  ); ++ i )
	{
	    const double & k = k_buf[ i ];
	    const double Rn = R_outer_integration( n, k );
	    p[ n ]->push_back( Rn );
	}		
    }

    std::cout << "Done." << std::endl;
    
    return;
}

double k_func::R_inner_integration
( int n, const double & r )
{
    if( r < nearly_0 )
	return 0.;
    else if( r > nearly_inf )
	switch( n )
	{
	case 1:
	    return 16. / 15.;
	case 2:
	    return -4. / 15.;
	default:
	    throw "R kernel error.";
	}
    else if( fabs( r - 1. ) < nearly_0 )
	switch( n )
	{
	case 1:
	    return 2. / 3.;
	case 2:
	    return 0.;
	default:
	    throw "R kernel error.";
	}
    else
	switch( n )
	{
	// This section looks messy:
	// I was using Mathematica "CForm" to get the
	// expressions...
	case 1:
	    return (4*r*(-3 + 11*pow(r,2) + 11*pow(r,4) 
			 - 3*pow(r,6))
		    - 3*pow(-1 + pow(r,2),4)
		*log(pow(-1 + r,2))
		    + 3*pow(-1 + pow(r,2),4)
		*log(pow(1 + r,2)))
		/(96.*pow(r,3));
	case 2:
	    return ((-1 + pow(r,2))
		    *(-4*r*(3 - 2*pow(r,2) + 3*pow(r,4))
		      - 3*pow(-1 + pow(r,2),2)
			*(1 + pow(r,2))*
		      log(pow(-1 + r,2))
			+ 3*pow(-1 + pow(r,2),2)
		      *(1 + pow(r,2))*log(pow(1 + r,2))))
		/(96.*pow(r,3));
	default:
	    throw  "R kernel error.";
	}
    return 0.;
}

double k_func::R_outer_integration
( int n, const double & k )
{
    integral intg;
    intg.clear(  );
    for( unsigned i = 0; i < k_buf.size(  ); ++ i )
    {
	const double r = k_buf[ i ] / k;
	const double kernel = PL_buf[ i ]
	    * R_inner_integration( n, r );
	intg.read( r, kernel );
    }

    static const double temp_coef = 0.025330295910584444;
    // 0.25 / \pi^2;
    return pow( k, 3 ) * temp_coef
	* PL_val( k ) * intg.result(  );
}

////////////////////////////////////////////////////////////
// Test output

void k_func::save_k_func( std::string file_name )
{
    std::ofstream fout( file_name.c_str(  ) );
    dvec * p[ 11 ] 
	= { &k_buf,   &R_buf_1, &R_buf_2,
	    &Q_buf_1, &Q_buf_2, &Q_buf_3,
	    &Q_buf_4, &Q_buf_5,	&Q_buf_6,
	    &Q_buf_7, &Q_buf_8 };

    for( unsigned i = 0; i < k_buf.size(  ); ++ i )
    {
	for( int j = 0; j < 11; ++ j )
	    fout << p[ j ]->at( i ) << ' ';
	fout << '\n';
    }
    fout << std::endl;

    return;
}

void k_func::load_k_func( std::string file_name )
{
    std::ifstream fin( file_name.c_str(  ) );
    if( !fin )
	throw "Unable to open k function file.";
	
    dvec * p[ 10 ] 
	= { &R_buf_1, &R_buf_2,
	    &Q_buf_1, &Q_buf_2, &Q_buf_3,
	    &Q_buf_4, &Q_buf_5,	&Q_buf_6,
	    &Q_buf_7, &Q_buf_8 };
    for( int i = 0; i < 10; ++ i )
	p[ i ]->clear(  );

    double temp( 0. );
    while( !fin.eof(  ) )
    {
	fin >> temp;
	for( int i = 0; i < 10; ++ i )
	{
	    fin >> temp;
	    p[ i ]->push_back( temp );
	}
    }
    for( int i = 1; i < 10; ++ i )
	p[ i ]->pop_back(  );

    std::cout << "k functions loaded from: "
	      << file_name << std::endl;
    return;
}
