#include "corr_func.h"
#include "lu_decomp.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

////////////////////////////////////////////////////////////
// Static variables

const double	corr_func::pi	    = 3.141592653589793;
const double	corr_func::max_y    = 100;
const int	corr_func::num_y    = 200;

////////////////////////////////////////////////////////////
// Constructor, destructor and initializer

corr_func::corr_func(  )
{

}

corr_func::~corr_func(  )
{
	
}
	
void corr_func::set_par( const corr_func_init & c_arg,
			 const q_func & qf )
{
    this->r_max	       = c_arg.r_max;
    this->r_min	       = c_arg.r_min;
    this->r_bin_num    = c_arg.r_bin_num;
    this->xi_file_name = c_arg.file_name;
    this->qf	       = ( q_func * )( &qf );

    std::cout << std::endl;
    return;
}

////////////////////////////////////////////////////////////
// Integration kernel


void corr_func::M( const double & r, const vec3 & y )
{
    vec3 q = y;
    q.z   += r;
    double q_norm = sqrt( q.x*q.x + q.y*q.y + q.z*q.z );
    const double qh[ 3 ]	// "q hat", unit vector
	= { q.x / q_norm, q.y / q_norm, q.z / q_norm };
    const double y_vec[ 3 ] = { y.x, y.y, y.z };
    
    q_func_vals qfv;
    qf->var_func( q_norm, qfv );
	
    const double xi_R = qfv.xi_L;

    ////////// Vectors and tensors //////////

    double A_temp[ 9 ];
    double A_inv[ 3 ][ 3 ], A_inv_temp[ 9 ];
    double A_det( 0. );
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    A_temp[lu.idx( i, j )]
		= delta_k( i, j )
		* ( qfv.X_11 + qfv.X_22 + 2.*qfv.X_13 )
		+ qh[ i ] * qh[ j ]
		* ( qfv.Y_11 + qfv.Y_22 + 2.*qfv.Y_13 );
    A_det = lu.lu_inverse( A_temp, A_inv_temp );
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    A_inv[ i ][ j ] = A_inv_temp[ lu.idx( i, j ) ];

    double A_10[ 3 ][ 3 ];
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    A_10[ i ][ j ]
		= delta_k( i, j ) * qfv.X_12_10
		+ qh[ i ] * qh[ j ] * qfv.Y_12_10
		+ delta_k( j, i ) * qfv.X_12_10
		+ qh[ j ] * qh[ i ] * qfv.Y_12_10;

    double U[ 3 ], U_20[ 3 ], U_11[ 3 ];
    for( int i = 0; i < 3; ++ i )
    {
	U[ i ] = ( qfv.U_1 + qfv.U_3 ) * qh[ i ];
	U_20[ i ] = qfv.U_2_20 * qh[ i ];
	U_11[ i ] = qfv.U_2_11 * qh[ i ];
    }

    double W[ 3 ][ 3 ][ 3 ], W_temp[ 3 ][ 3 ][ 3 ];
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    for( int k = 0; k < 3; ++ k )
		W_temp[ i ][ j ][ k ]
		    = qfv.V_112_1 * qh[ i ]
		      * delta_k( j, k )
		    + qfv.V_112_1 * qh[ j ]
		      * delta_k( k, i )
		    + qfv.V_112_3 * qh[ k ]
		      * delta_k( i, j )
		    + qfv.T_112
		      * qh[ i ] * qh[ j ] * qh[ k ];
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    for( int k = 0; k < 3; ++ k )
		W[ i ][ j ][ k ]
		    = W_temp[ i ][ j ][ k ]
		    + W_temp[ j ][ k ][ i ]
		    + W_temp[ k ][ i ][ j ];
	
    double g[ 3 ];
    for( int i = 0; i < 3; ++ i )
    {
	g[ i ] = 0.;
	for( int j = 0; j < 3; ++ j )
	    g[ i ] += A_inv[ i ][ j ] * y_vec[ j ];
    }

    double G[ 3 ][ 3 ];
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    G[ i ][ j ] = A_inv[ i ][ j ]
		- g[ i ] * g[ j ];

    double Gamma[ 3 ][ 3 ][ 3 ];
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    for( int k = 0; k < 3; ++ k )
		Gamma[ i ][ j ][ k ]
		    = A_inv[ j ][ k ] * g[ i ]
		    + A_inv[ k ][ i ] * g[ j ]
		    + A_inv[ i ][ j ] * g[ k ]
		    - g[ i ] * g[ j ] * g[ k ];

    ////////// Sum them up! //////////

    for( int i = 0; i < num_bias_comp; ++ i )
	bias_comp_inner[ i ] = 0.;
	
    double temp( 0. );
    bias_comp_inner[ 0 ] += 1.;

    bias_comp_inner[ 3 ] += xi_R;

    bias_comp_inner[ 5 ] += 0.5 * pow( xi_R, 2 );

    temp = 0.;
    for( int i = 0; i < 3; ++ i )
	temp += U[ i ] * g[ i ];
    bias_comp_inner[ 1 ] -= 2 * temp;
	
    temp = 0.;
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    temp += U[ i ] * U[ j ] * G[ i ][ j ];
    bias_comp_inner[ 2 ] -= temp;
    bias_comp_inner[ 3 ] -= temp;

    temp = 0.;
    for( int i = 0; i < 3; ++ i )
	temp += U[ i ] * g[ i ];
    bias_comp_inner[ 4 ] -= 2. * xi_R * temp;

    temp = 0.;
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    for( int k = 0; k < 3; ++ k )
		temp += W[ i ][ j ][ k ]
		    * Gamma[ i ][ j ][ k ];
    bias_comp_inner[ 0 ] += temp / 6.;

    temp = 0.;
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    temp += A_10[ i ][ j ] * G[ i ][ j ];
    bias_comp_inner[ 1 ] -= temp;

    temp = 0.;
    for( int i = 0; i < 3; ++ i )
	temp += U_20[ i ] * g[ i ];
    bias_comp_inner[ 2 ] -= temp;

    temp = 0.;
    for( int i = 0; i < 3; ++ i )
	temp += U_11[ i ] * g[ i ];
    bias_comp_inner[ 3 ] -= temp;

    temp = 0.;
    for( int i = 0; i < 3; ++ i )
	temp += y_vec[ i ] * g[ i ];
    static const double two_pi_cube = 248.05021344239853;
    // ( 2 \pi )^3
    const double gauss
	= exp( -0.5 * temp )
	/ sqrt( two_pi_cube * fabs( A_det ) );
    for( int i = 0; i < num_bias_comp; ++ i )
	bias_comp_inner[ i ] *= gauss;
	
    return;
}

void corr_func::M( const double & r,   const double & y,
		   const double & beta )
{
    vec3 y_vec;
    y_vec.x = y * sqrt( 1 - beta * beta );
    y_vec.y = 0;
    y_vec.z = y * beta;
    M( r, y_vec );
    return;
}

int corr_func::delta_k( const int & i, const int & j )
{
    return ( i == j ? 1 : 0 );
}

////////////////////////////////////////////////////////////
// Correlation function xi

void corr_func::xi( const double & r )
{
    double y( 0. ), beta( 0. );

    integral intg[ num_bias_comp ];
    for( int i = 0; i < num_bias_comp; ++ i )
	intg[ i ].clear(  );

    const double dy = max_y / num_y;
    while( y < max_y )
    {
	for( int k = 0; k < num_bias_comp; ++ k )
	    intg[ k ].gl_clear(  );
	for( int j = 0; j < intg[ 0 ].gl_num; ++ j )
	{
	    beta = intg[ 0 ].gl_xi( j );
	    M( r, y, beta );
	    for( int k = 0; k < num_bias_comp; ++ k )
		intg[ k ].gl_read
		    ( j, bias_comp_inner[ k ] );
	}
	for( int k = 0; k < num_bias_comp; ++ k )
	    intg[ k ].read
		( y, pow( y, 2 ) * intg[ k ].gl_result() );
	y += dy;
    }

    for( int k = 0; k < num_bias_comp; ++ k )
	bias_comp_outer[ k ]
	    = 2 * pi * intg[ k ].result(  );
    bias_comp_outer[ 0 ] -= 1.;
    
    return;
}

double corr_func::xi_L( const double & r )
{
    static q_func_vals qfv;
    qf->var_func( r, qfv );
	
    return qfv.xi_L;
}

void corr_func::get_xi(  )
{
    double r( 0. );
    std::cout << "Generate xi: ";
    std::cout.flush(  );

    std::vector<double> * p[ num_bias_comp ]
	={& b10b20, & b11b20, & b10b21,
	  & b12b20, & b11b21, & b10b22};

    pg.init( r_bin_num );
    for( int i = 0; i < r_bin_num; ++ i )
    {
	pg.show( i );
	r = ( r_max - r_min ) / ( r_bin_num - 1. )
	    * float( i ) + r_min;
	xi( r );
	r_buf.push_back( r );
	for( int j = 0; j < num_bias_comp; ++ j )
	    p[ j ]->push_back( bias_comp_outer[ j ] );
	xi_L_buf.push_back( xi_L( r ) );
    }

    std::cout << "Correlation function obtained. "
	      << std::endl;
    return;
}

////////////////////////////////////////////////////////////
// Output

void corr_func::output(  )
{
    std::ofstream fout( xi_file_name.c_str(  ) );
    std::vector<double> * p[ num_bias_comp + 2 ]
	={& r_buf, & xi_L_buf,
	  & b10b20, & b11b20, & b10b21,
	  & b12b20, & b11b21, & b10b22};
	
    for( unsigned i = 0; i < r_buf.size(  ); ++ i )
    {
	for( int j = 0; j < 8; ++ j )
	    fout << p[ j ]->at( i ) << ' ';
	fout << '\n';
    }
    fout.flush(  );

    return;
}


