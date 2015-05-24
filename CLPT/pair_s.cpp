#include "pair_s.h"
#include "lu_decomp.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

////////////////////////////////////////////////////////////
// Static variables

const double	pair_s::pi	 = 3.141592653589793;
const double	pair_s::max_y    = 100;
const int	pair_s::num_y    = 200;

////////////////////////////////////////////////////////////
// Constructor, destructor and initializer

pair_s::pair_s(  )
{
	
}

pair_s::~pair_s(  )
{
	
}

void pair_s::set_par( const corr_func_init & s_arg,
                      const q_func & qf )
{
    this->r_max		= s_arg.r_max;
    this->r_min		= s_arg.r_min;
    this->r_bin_num	= s_arg.r_bin_num;
    this->s12_file_name = s_arg.file_name;
    this->qf		= ( q_func * )( & qf );

    std::cout << std::endl;
    return;
}

////////////////////////////////////////////////////////////
// Integration kernel

void pair_s::M2( const double & r, const vec3 & y )
{
    // Direction index; for testing currently.
    static const double rh[ 3 ]  = { 0, 0, 1 };
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
	
    // for( int i = 0; i < num_bias_comp; ++ i )
    // {
    // 	bias_comp_inner_p[ i ] = 0.;
    // 	bias_comp_inner_v[ i ] = 0.;
    // }

    double temp_p( 0. ), temp_v( 0. );
    // parallel and perpandicular (vertical)
    // components of the result.
    double temp[ 3 ][ 3 ];
    // I would sacrifise a little efficiency
    // so that the code does not look too bad...

    // b1^0 b2^0 -> bias_comp[ 0 ] -> b0
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	    temp[ n ][ m ] = A_ddot[ n ][ m ];
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	    for( int i = 0; i < 3; ++ i )
		temp[ n ][ m ]
		    -= g[ i ] * W_ddot[ i ][ n ][ m ];
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	    for( int i = 0; i < 3; ++ i )
		for( int j = 0; j < 3; ++ j )
		    temp[ n ][ m ]
			-= A_dot[ i ][ n ] 
			 * A_dot[ j ][ m ] * G[ i ][ j ];
    temp_p = 0.;
    temp_v = 0.;
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	{
	    temp_p += temp[ n ][ m ]
		* rh[ n ] * rh[ m ];
	    temp_v += temp[ n ][ m ]
		* ( lh1[ n ] * lh1[ m ]
		    + lh2[ n ] * lh2[ m ]);
	}
    bias_comp_inner_p[ 0 ] = temp_p;
    bias_comp_inner_v[ 0 ] = temp_v;

    // b1^1 b2^0 -> bias_comp[ 1 ] -> b11
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	    temp[ n ][ m ] = 2. * A_ddot_10[ n ][ m ];
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	    for( int i = 0; i < 3; ++ i )
	    {
		temp[ n ][ m ] -= 2. *
		    ( A_dot[ i ][ n ]
			* g[ i ] * U_dot[ m ]
		    + A_dot[ i ][ m ]
			* g[ i ] * U_dot[ n ] );
		temp[ n ][ m ] -= 2. *
		    U[ i ] * g[ i ] * A_ddot[ n ][ m ];
	    }
    temp_p = 0.;
    temp_v = 0.;
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	{
	    temp_p += temp[ n ][ m ]
		* rh[ n ] * rh[ m ];
	    temp_v += temp[ n ][ m ]
		* ( lh1[ n ] * lh1[ m ]
		    + lh2[ n ] * lh2[ m ]);
	}
    bias_comp_inner_p[ 1 ] = temp_p;
    bias_comp_inner_v[ 1 ] = temp_v;

    // b1^0 b2^1 -> bias_comp[ 2 ] -> b21
    temp_p = 0.;
    temp_v = 0.;
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	    temp[ n ][ m ] = 2. * U_dot[ n ] * U_dot[ m ];
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	{
	    temp_p += temp[ n ][ m ]
		* rh[ n ] * rh[ m ];
	    temp_v += temp[ n ][ m ]
		* ( lh1[ n ] * lh1[ m ]
		    + lh2[ n ] * lh2[ m ]);
	}
    bias_comp_inner_p[ 2 ] = temp_p;
    bias_comp_inner_v[ 2 ] = temp_v;

    // b1^2 b2^0 -> bias_comp[ 3 ] -> b12
    temp_p = 0.;
    temp_v = 0.;
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	    temp[ n ][ m ]
		= 2. * U_dot[ n ] * U_dot[ m ]
		+ xi_R * A_ddot[ n ][ m ];
    for( int n = 0; n < 3; ++ n )
	for( int m = 0; m < 3; ++ m )
	{
	    temp_p += temp[ n ][ m ]
		* rh[ n ] * rh[ m ];
	    temp_v += temp[ n ][ m ]
		* ( lh1[ n ] * lh1[ m ]
		    + lh2[ n ] * lh2[ m ]);
	}
    bias_comp_inner_p[ 3 ] = temp_p;
    bias_comp_inner_v[ 3 ] = temp_v;

    // Gaussian-like factor
    temp_p = 0.;
    for( int i = 0; i < 3; ++ i )
	temp_p += y_vec[ i ] * g[ i ];

    static const double two_pi_cube = 248.05021344239853;
    // ( 2 \pi )^3
    const double gauss
	= exp( -0.5 * temp_p )
	/ sqrt( two_pi_cube * fabs( A_det ) );
    for( int i = 0; i < num_bias_comp; ++ i )
    {
	bias_comp_inner_p[ i ] *= gauss;
	bias_comp_inner_v[ i ] *= gauss;
    }
	
    return;
}

void pair_s::M2( const double & r, const double & y,
                 const double & beta )
{
    vec3 y_vec;
    y_vec.x = y * sqrt( 1 - beta * beta );
    y_vec.y = 0;
    y_vec.z = y * beta;
    M2( r, y_vec );
    return;
}

int pair_s::delta_k( const int & i, const int & j )
{
    return ( i == j ? 1 : 0 );
}

////////////////////////////////////////////////////////////
// Correlation function v12

void pair_s::s12( const double & r )
{
    double y( 0. ), beta( 0. );

    integral intg_p[ num_bias_comp ];
    integral intg_v[ num_bias_comp ];
    for( int i = 0; i < num_bias_comp; ++ i )
    {
	intg_p[ i ].clear(  );
	intg_v[ i ].clear(  );
    }

    const double dy = max_y / num_y;
    while( y < max_y )	
    {
	for( int k = 0; k < num_bias_comp; ++ k )
	{
	    intg_p[ k ].gl_clear(  );
	    intg_v[ k ].gl_clear(  );
	}
	for( int j = 0; j < intg_p[ 0 ].gl_num; ++ j )
	{
	    beta = intg_p[ 0 ].gl_xi( j );
	    M2( r, y, beta );
	    for( int k = 0; k < num_bias_comp; ++ k )
	    {
		intg_p[ k ].gl_read
		    ( j, bias_comp_inner_p[ k ] );
		intg_v[ k ].gl_read
		    ( j, bias_comp_inner_v[ k ] );
	    }
	}
	for( int k = 0; k < num_bias_comp; ++ k )
	{
	    intg_p[ k ].read( y, pow( y, 2 )
		              * intg_p[ k ].gl_result(  ) );
	    intg_v[ k ].read( y, pow( y, 2 )
 		              * intg_v[ k ].gl_result(  ) );
	}
	y += dy;
    }

    for( int k = 0; k < num_bias_comp; ++ k )
    {
	bias_comp_outer_p[ k ]
	    = 2 * pi * intg_p[ k ].result(  );
	bias_comp_outer_v[ k ]
	    = 2 * pi * intg_v[ k ].result(  );
    }
    return;
}

void pair_s::get_s12(  )
{
    double r( 0. );
    std::cout << "Generate sigma12: ";
    std::cout.flush(  );

    std::vector<double> * p_p[ num_bias_comp ]
	={ & b0p, & b11p, & b21p, & b12p };
    std::vector<double> * p_v[ num_bias_comp ]
	={ & b0v, & b11v, & b21v, & b12v };

    pg.init( r_bin_num );
    const double dr = ( r_max - r_min )
	/ ( r_bin_num - 1. );
	
    for( int i = 0; i < r_bin_num; ++ i )
    {
	pg.show( i );
	r = dr * i + r_min;
	s12( r );
	r_buf.push_back( r );
	for( int j = 0; j < num_bias_comp; ++ j )
	{
	    p_p[ j ]->push_back( bias_comp_outer_p[ j ] );
	    p_v[ j ]->push_back( bias_comp_outer_v[ j ] );
	}
    }

    std::cout << "Pairwise velocity dispersion obtained. "
	      << std::endl;
    return;
}

////////////////////////////////////////////////////////////
// Output

void pair_s::output(  )
{
    std::string temp_file_name;
    temp_file_name = s12_file_name + "_p";
    std::ofstream fout_p( temp_file_name.c_str(  ) );
    temp_file_name = s12_file_name + "_v";
    std::ofstream fout_v( temp_file_name.c_str(  ) );
	
    std::vector<double> * p_p[ num_bias_comp ]
	={ & b0p, & b11p, & b21p, & b12p };
    std::vector<double> * p_v[ num_bias_comp ]
	={ & b0v, & b11v, & b21v, & b12v };
	
    for( unsigned i = 0; i < r_buf.size(  ); ++ i )
    {
	fout_p << r_buf[ i ] << ' ';
	fout_v << r_buf[ i ] << ' ';
	for( int j = 0; j < num_bias_comp; ++ j )
	{
	    fout_p << p_p[ j ]->at( i ) << ' ';
	    fout_v << p_v[ j ]->at( i ) << ' ';			
	}
	fout_p << '\n';
	fout_v << '\n';
    }
	
    return;
}


