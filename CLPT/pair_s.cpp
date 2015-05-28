#include "pair_s.h"
#include "lu_decomp.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

////////////////////////////////////////////////////////////
// Constructor, destructor and initializer

pair_s::pair_s(  )
{
    num_bias_comp = 8;
}

////////////////////////////////////////////////////////////
// Integration kernel

void pair_s::kernel( const double & r, const vec3 & y,
                     double * bias_comp )
{
    // The parallel component comes first
    double * bias_comp_p = bias_comp;
    double * bias_comp_v = bias_comp + 4;
    
    // Direction index; for testing currently.
    static const double rh [ 3 ] = { 0, 0, 1 };
    static const double lh2[ 3 ] = { 0, 1, 0 };
    static const double lh1[ 3 ] = { 1, 0, 0 };

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
	    A_inv[ i ][ j ] = A_inv_temp[lu.idx( i, j )];

    double U[ 3 ];
    for( int i = 0; i < 3; ++ i )
	U[ i ] = ( qfv.U_1 + qfv.U_3 ) * qh[ i ];

    double g[ 3 ];
    for( int i = 0; i < 3; ++ i )
	g[ i ] = 0.;
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    g[ i ] += A_inv[ i ][ j ] * y_vec[ j ];

    double G[ 3 ][ 3 ];
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    G[ i ][ j ] = A_inv[ i ][ j ]
		- g[ i ] * g[ j ];

    double A_dot[ 3 ][ 3 ];
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    A_dot[ i ][ j ] = delta_k( i, j )
		* ( qfv.X_11 + 2. * qfv.X_22
		    + 4. * qfv.X_13 )
		+ qh[ i ] * qh[ j ]
		* ( qfv.Y_11 + 2. * qfv.Y_22
		    + 4. * qfv.Y_13 );

    double A_ddot[ 3 ][ 3 ];
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    A_ddot[ i ][ j ] = delta_k( i, j )
		* ( qfv.X_11 + 4. * qfv.X_22
		    + 6. * qfv.X_13 )
		+ qh[ i ] * qh[ j ]
		* ( qfv.Y_11 + 4. * qfv.Y_22
		    + 6. * qfv.Y_13 );

    double A_ddot_10[ 3 ][ 3 ];
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    A_ddot_10[ i ][ j ]	= 4. *
		( qfv.X_12_10 * delta_k( i, j )
		    + qfv.Y_12_10 * qh[ i ] * qh[ j ] );

    double U_dot[ 3 ];
    for( int i = 0; i < 3; ++ i )
	U_dot[ i ] = qh[ i ] * ( qfv.U_1 + 3. * qfv.U_3 );

    double W_temp[ 3 ][ 3 ][ 3 ];
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

    double W_ddot[ 3 ][ 3 ][ 3 ];
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    for( int k = 0; k < 3; ++ k )
		W_ddot[ i ][ j ][ k ]
		    = 2. * W_temp[ i ][ j ][ k ]
		    + 2. * W_temp[ i ][ k ][ j ]
		    + W_temp[ k ][ j ][ i ];


    ////////// Sum them up! //////////
    double sum_p( 0. ), sum_v( 0. );
    // parallel and perpandicular (vertical)
    // components of the result.
    double sum[ 3 ][ 3 ];
    // I would sacrifise a little efficiency
    // so that the code does not look too bad...

    // b1^0 b2^0 -> bias_comp[ 0 ] -> b0
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	    sum[ n ][ m ] = A_ddot[ n ][ m ];
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	    for( int i = 0; i < 3; ++ i )
		sum[ n ][ m ]
		    -= g[ i ] * W_ddot[ i ][ n ][ m ];
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	    for( int i = 0; i < 3; ++ i )
		for( int j = 0; j < 3; ++ j )
		    sum[ n ][ m ]
			-= A_dot[ i ][ n ]
			 * A_dot[ j ][ m ] * G[ i ][ j ];
    sum_p = 0.;
    sum_v = 0.;
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	{
	    sum_p += sum[ n ][ m ]
		* rh[ n ] * rh[ m ];
	    sum_v += sum[ n ][ m ]
		* ( lh1[ n ] * lh1[ m ]
		    + lh2[ n ] * lh2[ m ]);
	}
    bias_comp_p[ 0 ] = sum_p;
    bias_comp_v[ 0 ] = sum_v;

    // b1^1 b2^0 -> bias_comp[ 1 ] -> b11
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	    sum[ n ][ m ] = 2. * A_ddot_10[ n ][ m ];
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	    for( int i = 0; i < 3; ++ i )
	    {
		sum[ n ][ m ] -= 2. *
		    ( A_dot[ i ][ n ]
			* g[ i ] * U_dot[ m ]
		    + A_dot[ i ][ m ]
			* g[ i ] * U_dot[ n ] );
		sum[ n ][ m ] -= 2. *
		    U[ i ] * g[ i ] * A_ddot[ n ][ m ];
	    }
    sum_p = 0.;
    sum_v = 0.;
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	{
	    sum_p += sum[ n ][ m ]
		* rh[ n ] * rh[ m ];
	    sum_v += sum[ n ][ m ]
		* ( lh1[ n ] * lh1[ m ]
		    + lh2[ n ] * lh2[ m ]);
	}
    bias_comp_p[ 1 ] = sum_p;
    bias_comp_v[ 1 ] = sum_v;

    // b1^0 b2^1 -> bias_comp[ 2 ] -> b21
    sum_p = 0.;
    sum_v = 0.;
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	    sum[ n ][ m ] = 2. * U_dot[ n ] * U_dot[ m ];
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	{
	    sum_p += sum[ n ][ m ]
		* rh[ n ] * rh[ m ];
	    sum_v += sum[ n ][ m ]
		* ( lh1[ n ] * lh1[ m ]
		    + lh2[ n ] * lh2[ m ]);
	}
    bias_comp_p[ 2 ] = sum_p;
    bias_comp_v[ 2 ] = sum_v;

    // b1^2 b2^0 -> bias_comp[ 3 ] -> b12
    sum_p = 0.;
    sum_v = 0.;
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	    sum[ n ][ m ]
		= 2. * U_dot[ n ] * U_dot[ m ]
		+ xi_R * A_ddot[ n ][ m ];
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	{
	    sum_p += sum[ n ][ m ]
		* rh[ n ] * rh[ m ];
	    sum_v += sum[ n ][ m ]
		* ( lh1[ n ] * lh1[ m ]
		    + lh2[ n ] * lh2[ m ]);
	}
    bias_comp_p[ 3 ] = sum_p;
    bias_comp_v[ 3 ] = sum_v;

    // Gaussian-like factor
    double gauss_exp( 0. );
    for( int i = 0; i < 3; ++ i )
	gauss_exp += y_vec[ i ] * g[ i ];
    static const double two_pi_cube = 248.05021344239853;
    // ( 2 \pi )^3
    const double gauss = exp( -0.5 * gauss_exp )
	/ sqrt( two_pi_cube * fabs( A_det ) );
    for( unsigned i = 0; i < num_bias_comp; ++ i )
	bias_comp[ i ] *= gauss;

    return;
}

void pair_s::post_proc(  )
{
    std::cout << "Pairwise sigma^2";
    return;
}
