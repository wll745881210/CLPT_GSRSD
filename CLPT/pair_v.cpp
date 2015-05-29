#include "pair_v.h"

#include <iostream>
#include <cmath>

pair_v::pair_v(  )
{
    num_bias_comp = 5;
}

void pair_v::kernel( const double & r, const vec3 & y,
                     double * bias_comp )
{
    // Direction index; for testing currently.
    const double rh[ 3 ] = { 0, 0, 1 };

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
		* ( qfv.X_11 + qfv.X_22 + 2 * qfv.X_13 )
		+ qh[ i ] * qh[ j ]
		* ( qfv.Y_11 + qfv.Y_22 + 2 * qfv.Y_13 );
    A_det = lu.lu_inverse( A_temp, A_inv_temp );
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    A_inv[ i ][ j ] = A_inv_temp[lu.idx( i, j )];

    double U[ 3 ];
    for( int i = 0; i < 3; ++ i )
	U[ i ] = ( qfv.U_1 + qfv.U_3 ) * qh[ i ];
	
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
	
    double A_dot[ 3 ][ 3 ];
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    A_dot[ i ][ j ] = delta_k( i, j )
		* ( qfv.X_11 + 2. * qfv.X_22
		    + 4. * qfv.X_13 )
		+ qh[ i ] * qh[ j ]
		* ( qfv.Y_11 + 2. * qfv.Y_22
		    + 4. * qfv.Y_13 );

    double U_dot[ 3 ], U_20_dot[ 3 ], U_11_dot[ 3 ];
    for( int i = 0; i < 3; ++ i )
    {
	U_dot[ i ] = qh[ i ] * ( qfv.U_1 + 3. * qfv.U_3 );
	U_20_dot[ i ] = 2. * qh[ i ] * qfv.U_2_20;
	U_11_dot[ i ] = 2. * qh[ i ] * qfv.U_2_11;
    }

    double W_dot[ 3 ][ 3 ][ 3 ];
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    for( int k = 0; k < 3; ++ k )
		W_dot[ i ][ j ][ k ]
		    = ( 3. * qfv.V_112_1 + qfv.V_112_3 )
		    * ( qh[ i ] * delta_k( j, k )
			+ qh[ j ] * delta_k( k, i ) )
		    + 2. * ( qfv.V_112_1 + qfv.V_112_3 )
		    * qh[ k ] * delta_k( i, j )
		    + 4. * qfv.T_112
		    * qh[ i ] * qh[ j ] * qh[ k ];

    double A_dot_10[ 3 ][ 3 ];
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    A_dot_10[ i ][ j ]
		= 3. *   // Should be 3, not 2 or 4 ?
		( qfv.X_12_10 * delta_k( i, j )
		+ qfv.Y_12_10 * qh[ i ] * qh[ j ] );
	

    ////////// Sum them up! //////////
	
    for( unsigned i = 0; i < num_bias_comp; ++ i )
	bias_comp[ i ] = 0.;

    double sum( 0. );
	
    sum = 0.;
    for( int i = 0; i < 3; ++ i )
	for( int n = 0; n < 3; ++ n )
	    sum += g[ i ] * rh[ n ]
		* A_dot[ i ][ n ];
    bias_comp[ 0 ] -= sum;

    sum = 0.;
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    for( int n = 0; n < 3; ++ n )
		sum += rh[ n ] * G[ i ][ j ]
		    * W_dot[ i ][ j ][ n ];
    bias_comp[ 0 ] -= 0.5 * sum;

    sum = 0.;
    for( int n = 0; n < 3; ++ n )
	sum += U_dot[ n ] * rh[ n ];
    bias_comp[ 1 ] += 2. * sum;
	
    sum = 0.;
    for( int i = 0; i < 3; ++ i )
	for( int n = 0; n < 3; ++ n )
	    sum += g[ i ] * A_dot_10[ i ][ n ]
		* rh[ n ];
    bias_comp[ 1 ] -= 2. * sum;

    sum = 0.;
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    for( int n = 0; n < 3; ++ n )
		sum += rh[ n ] * G[ i ][ j ]
		    * U[ i ] * A_dot[ j ][ n ];
    bias_comp[ 1 ] -= 2. * sum;


    sum = 0.;
    for( int n = 0; n < 3; ++ n )
	sum += U_20_dot[ n ] * rh[ n ];
    bias_comp[ 2 ] += sum;

    sum = 0.;
    for( int i = 0; i < 3; ++ i )
	for( int n = 0; n < 3; ++ n )
	    sum += g[ i ] * U[ i ]
		* U_dot[ n ] * rh[ n ];
    bias_comp[ 2 ] -= 2. * sum;
    bias_comp[ 3 ] -= 2. * sum;

    sum = 0.;
    for( int i = 0; i < 3; ++ i )
	for( int n = 0; n < 3; ++ n )
	    sum += g[ i ] * rh[ n ]
		* A_dot[ i ][ n ];
    bias_comp[ 3 ] -= xi_R * sum;

    sum = 0.;
    for( int n = 0; n < 3; ++ n )
	sum += U_11_dot[ n ] * rh[ n ];
    bias_comp[ 3 ] += sum;


    sum = 0.;
    for( int n = 0; n < 3; ++ n )
	sum += rh[ n ] * U_dot[ n ];
    bias_comp[ 4 ] += 2. * xi_R * sum;
	
    sum = 0.;
    for( int i = 0; i < 3; ++ i )
	sum += y[ i ] * g[ i ];

    static const double two_pi_cube = 248.05021344239853;
    // ( 2 \pi )^3
    const double gauss
	= exp( -0.5 * sum )
	/ sqrt( two_pi_cube * fabs( A_det ) );
    for( unsigned i = 0; i < num_bias_comp; ++ i )
	bias_comp[ i ] *= gauss;
	
    return;
}

void pair_v::post_proc(  )
{
    std::vector<double> * v12_L_vec
	= new std::vector<double> ( rvec.size(  ) );
#pragma omp parallel for
    for( unsigned i = 0; i < rvec.size(  ); ++ i )
	v12_L_vec->at( i ) = this->v12_L( rvec[ i ] );

    corr_res.insert( corr_res.begin(  ), v12_L_vec );
    std::cout << "Pairwise velocity";
    return;
}


double pair_v::v12_L( const double & r )
{
    integral intg;
    intg.clear(  );
    k_func * kf = k_func::get_instance(  );
	
    const std::vector<double> & kv = kf->kvec(  );

    for( unsigned i = 0; i < kv.size(  ); ++ i )
    {
	const double & k = kv[ i ];
	const double  jx = k * r;
	const double kernel = k * kf->PL_val( k )
	    * ( sin( jx ) - jx * cos( jx ) )
	    / pow( jx, 2 );
	intg.read( k, kernel );
    }
    return -1. / pow( pi, 2 ) * intg.result(  );
}
