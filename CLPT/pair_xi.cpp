#include "pair_xi.h"

#include <iostream>
#include <cmath>

void pair_xi::kernel( const double & r, const vec3 & y,
                      double * bias_comp )
{
    vec3 q = y;
    q.z   += r;
    const vec3   qhat = q.vhat(  ); 	// unit vector
    const double qh[] = { qhat[ 0 ], qhat[ 1 ], qhat[ 2 ] };

    q_func_vals qfv;
    qf->var_func( q.norm(  ), qfv );
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
	    g[ i ] += A_inv[ i ][ j ] * y[ j ];
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

    for( unsigned i = 0; i < num_bias_comp; ++ i )
	bias_comp[ i ] = 0.;
	
    double sum( 0. );
    bias_comp[ 0 ] += 1.;

    bias_comp[ 3 ] += xi_R;

    bias_comp[ 5 ] += 0.5 * pow( xi_R, 2 );

    sum = 0.;
    for( int i = 0; i < 3; ++ i )
	sum += U[ i ] * g[ i ];
    bias_comp[ 1 ] -= 2 * sum;
	
    sum = 0.;
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    sum += U[ i ] * U[ j ] * G[ i ][ j ];
    bias_comp[ 2 ] -= sum;
    bias_comp[ 3 ] -= sum;

    sum = 0.;
    for( int i = 0; i < 3; ++ i )
	sum += U[ i ] * g[ i ];
    bias_comp[ 4 ] -= 2. * xi_R * sum;

    sum = 0.;
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    for( int k = 0; k < 3; ++ k )
		sum += W[ i ][ j ][ k ]
		    * Gamma[ i ][ j ][ k ];
    bias_comp[ 0 ] += sum / 6.;

    sum = 0.;
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    sum += A_10[ i ][ j ] * G[ i ][ j ];
    bias_comp[ 1 ] -= sum;

    sum = 0.;
    for( int i = 0; i < 3; ++ i )
	sum += U_20[ i ] * g[ i ];
    bias_comp[ 2 ] -= sum;

    sum = 0.;
    for( int i = 0; i < 3; ++ i )
	sum += U_11[ i ] * g[ i ];
    bias_comp[ 3 ] -= sum;

    sum = 0.;
    for( int i = 0; i < 3; ++ i )
	sum += y[ i ] * g[ i ];
    static const double two_pi_cube = 248.05021344239853;
    // ( 2 \pi )^3
    const double gauss = exp( -0.5 * sum )
	/ sqrt( two_pi_cube * fabs( A_det ) );
    for( unsigned i = 0; i < num_bias_comp; ++ i )
	bias_comp[ i ] *= gauss;
	
    return;
}

void pair_xi::post_proc(  )
{
    std::vector<double> * xi_L = new std::vector<double>;
    q_func_vals qval;
    
    for( unsigned i = 0; i < rvec.size(  ); ++ i )
    {
	corr_res[ 0 ]->at( i ) -= 1.;
	qf->var_func( rvec[ i ], qval );
	xi_L->push_back( qval.xi_L );
    }

    corr_res.insert( corr_res.begin(  ), xi_L );
    std::cout << "Correlation xi";
    return;
}

