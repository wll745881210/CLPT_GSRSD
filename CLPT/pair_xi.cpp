#include "pair_xi.h"

#include <iostream>
#include <cmath>

void pair_xi::kernel( const double & r, const vec3 & y )
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
