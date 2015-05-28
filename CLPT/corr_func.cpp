#include "corr_func.h"
#include "lu_decomp.h"
#include "save_load.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

////////////////////////////////////////////////////////////
// Static variables

const double corr_func::pi = 3.141592653589793;

double    corr_func::max_y    = 100;
int       corr_func::num_y    = 200;
q_func *  corr_func::qf       = NULL;

////////////////////////////////////////////////////////////
// Constructor, destructor and initializer

corr_func::corr_func(  )
{

}

corr_func::~corr_func(  )
{
	
}
	
void corr_func::set_par( input & args )
{
    double r_min( 0. ), r_max( 0. );
    int r_bin( 0 );
    args.find_key( "r_max",     r_max, 130 );
    args.find_key( "r_min",     r_max, 1   );
    args.find_key( "r_bin_num", r_bin, 80  );
    const double dr = ( r_max - r_min ) / ( r_bin - 1 );
    for( int i = 0; i < r_bin; ++ i  )
	rvec.push_back( r_min + i * dr );

    qf = q_func::get_instance(  );

    args.find_key

    return;
}

void corr_func::initialize(  )
{
    static const unsigned num_bias_comp = 6;
    bias_comp = new double[ num_bias_comp ];
    for( unsigned i = 0; i < num_bias_comp; ++ i )
    {
	std::vector<double> * p
	    = new std::vector<double> ( rvec.size(  ) );
	corr_res.push_back( p );
    }
    return;
}

////////////////////////////////////////////////////////////
// Integration kernel itself is virtual



void corr_func::ker_wrap
( const double & r, const double & y, const double & beta )
{
    vec3 y_vec;
    y_vec.x = y * sqrt( 1 - pow( beta, 2 ) );
    y_vec.y = 0;
    y_vec.z = y * beta;
    kernel( r, y_vec );
    return;
}

int corr_func::delta_k( const int & i, const int & j )
{
    return ( i == j ? 1 : 0 );
}

////////////////////////////////////////////////////////////
// Correlation function xi

void corr_func::calc_corr( const double & r )
{
    integral intg[ num_bias_comp ];
    for( int i = 0; i < num_bias_comp; ++ i )
	intg[ i ].clear(  );

    const double dy = max_y / num_y;
    for( double y = 0.; y < max_y; y += dy )
    {
	for( int k = 0; k < num_bias_comp; ++ k )
	    intg[ k ].gl_clear(  );
	for( int j = 0; j < intg[ 0 ].gl_num; ++ j )
	{
	    const double beta = intg[ 0 ].gl_xi( j );
	    M( r, y, beta );
	    for( int k = 0; k < num_bias_comp; ++ k )
		intg[ k ].gl_read
		    ( j, bias_comp_inner[ k ] );
	}
	for( int k = 0; k < num_bias_comp; ++ k )
	    intg[ k ].read
		( y, pow( y, 2 ) * intg[ k ].gl_result() );
    }

    for( int k = 0; k < num_bias_comp; ++ k )
	bias_comp_outer[ k ]
	    = 2 * pi * intg[ k ].result(  );
    bias_comp_outer[ 0 ] -= 1.;
    
    return;
}

void corr_func::get_corr(  )
{
    std::cout << "Generate xi: ";
    std::cout.flush(  );

    std::vector<double> * p[ num_bias_comp ]
	={& b10b20, & b11b20, & b10b21,
	  & b12b20, & b11b21, & b10b22};

    for( int i = 0; i < r_bin_num; ++ i )
	calc_corr( i );

    std::cout << "Correlation function obtained. "
	      << std::endl;
    return;
}

////////////////////////////////////////////////////////////
// Output

void corr_func::output( const std::string file_path )
{
    save_load s( file_path );
    for(  )
    return;
}


