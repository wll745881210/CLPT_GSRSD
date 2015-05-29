#include "corr_func.h"
#include "lu_decomp.h"
#include "save_load.h"
#include "q_depend_funcs.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

////////////////////////////////////////////////////////////
// For the vector

double vec3::norm2(  )
{
    return pow( x, 2 ) + pow( y, 2 ) + pow( z, 2 );
}

double vec3::norm(  )
{
    return sqrt( this->norm2(  ) );
}

vec3 vec3::vhat(  )
{
    vec3 res = ( * this );
    const double norm_val = this->norm(  );
    for( int i = 0; i < 3; ++ i )
	res[ i ] /= norm_val;
    return res;
}

double & vec3::operator[] ( const unsigned & i )
{
    if( i == 0 )
	return this->x;
    else if( i == 1 )
	return this->y;
    else
	return this->z;
}

const double & vec3::operator[] ( const unsigned & i ) const
{
    if( i == 0 )
    	return this->x;
    else if( i == 1 )
    	return this->y;
    else
    	return this->z;
}

////////////////////////////////////////////////////////////
// Static variables

const double corr_func::pi( 3.141592653589793 );
double    corr_func::y_max( 25 );
int       corr_func::y_bin( 50 );
q_func *  corr_func::qf    = q_func::get_instance(  );

std::vector<double> corr_func::rvec;

////////////////////////////////////////////////////////////
// Constructor, destructor and initializer

corr_func::corr_func(  )
{
    num_bias_comp = 6;
}

corr_func::~corr_func(  )
{
    for( unsigned i = 0; i < corr_res.size(  ); ++ i )
	delete ( corr_res[ i ] );
    return;
}
	
void corr_func::set_par( input & args )
{
    double r_min( 0. ), r_max( 0. );
    int r_bin( 0 );
    args.find_key( "r_max",     r_max, 130 );
    args.find_key( "r_min",     r_min, 1   );
    args.find_key( "r_bin_num", r_bin, 80  );
    const double dr = ( r_max - r_min ) / ( r_bin - 1 );
    for( int i = 0; i < r_bin; ++ i  )
	rvec.push_back( r_min + i * dr );

    args.find_key( "y_max",      y_max, 50 );
    args.find_key( "y_bin_num",  y_bin, 100 );
 
    return;
}

void corr_func::initialize(  )
{
    for( unsigned i = 0; i < num_bias_comp; ++ i )
    {
	std::vector<double> * p
	    = new std::vector<double> ( rvec.size(  ) );
	corr_res.push_back( p );
    }
    return;
}

////////////////////////////////////////////////////////////
// Integration kernel is purely virtual; its wrapper is not

void corr_func::ker_wrap
( const double & r, const double & y, const double & beta,
  double * bias_comp )
{
    vec3 y_vec;
    y_vec.x = y * sqrt( 1 - pow( beta, 2 ) );
    y_vec.y = 0;
    y_vec.z = y * beta;
    return kernel( r, y_vec, bias_comp );
}

int corr_func::delta_k( const int & i, const int & j )
{
    return ( i == j ? 1 : 0 );
}

////////////////////////////////////////////////////////////
// Correlation function xi

void corr_func::calc_corr( const unsigned & i )
{
    const double & r = rvec[ i ];
    
    integral * intg = new integral[ num_bias_comp ];
    double bias_comp[ num_bias_comp ];
    const double dy = y_max / y_bin;
    for( double y = 0.; y < y_max; y += dy )
    {
	for( unsigned j = 0; j < num_bias_comp; ++ j )
	    intg[ j ].gl_clear(  );
	for( int l = 0; l < intg[ 0 ].gl_num; ++ l )
	{
	    const double beta = intg[ 0 ].gl_xi( l );
	    ker_wrap( r, y, beta, bias_comp );
	    for( unsigned j = 0; j < num_bias_comp; ++ j )
		intg[ j ].gl_read( l, bias_comp[ j ] );
	}
	for( unsigned j = 0; j < num_bias_comp; ++ j )
	{
	    const double inner = intg[ j ].gl_result(  );
	    intg[ j ].read( y, pow( y, 2 ) * inner );
	}
    }

    for( unsigned j = 0; j < num_bias_comp; ++ j )
	corr_res[ j ]->at( i )
	    = 2 * pi * intg[ j ].result(  );
    return;
}

void corr_func::post_proc(  )
{
    return;
}

void corr_func::get_corr(  )
{
    std::cout << "Generating correlation... " << std::flush;
#pragma omp parallel for schedule( dynamic, 10 )
    for( unsigned i = 0; i < rvec.size(  ); ++ i )
	calc_corr( i );

    post_proc(  );
    std::cout << " function obtained. "
	      << std::endl;
    return;
}

////////////////////////////////////////////////////////////
// Output

void corr_func::output( const std::string file_path )
{
    save_load s( file_path );
    s.add_vec( rvec );
    for( unsigned i = 0; i < corr_res.size(  ); ++ i )
	s.add_vec( * corr_res[ i ] );
    s.write(  );
    
    return;
}


