#include "pair_v.h"
#include "lu_decomp.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

////////////////////////////////////////////////////////////
// Static variables

const double	pair_v::pi	 = 3.141592653589793;
const double	pair_v::max_y    = 100;
const int	pair_v::num_y    = 200;

////////////////////////////////////////////////////////////
// Constructor, destructor and initializer

pair_v::pair_v(  )
{
	
}

pair_v::~pair_v(  )
{
	
}

void pair_v::set_par( const corr_func_init & v_arg,
                      const q_func & qf )
{
    this->r_max		= v_arg.r_max;
    this->r_min		= v_arg.r_min;
    this->r_bin_num	= v_arg.r_bin_num;
    this->v12_file_name = v_arg.file_name;
    this->qf		= ( q_func * )( & qf );

    std::cout << std::endl;
    return;
}

////////////////////////////////////////////////////////////
// Integration kernel

void pair_v::M1( const double & r, const vec3 & y )
{
    // Direction index; for testing currently.
    const double rh[ 3 ] = { 0, 0, 1 };

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
	    g[ i ] += A_inv[ i ][ j ] * y_vec[ j ];
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
	
    for( int i = 0; i < num_bias_comp; ++ i )
	bias_comp_inner[ i ] = 0.;

    double temp( 0. );
	
    temp = 0.;
    for( int i = 0; i < 3; ++ i )
	for( int n = 0; n < 3; ++ n )
	    temp += g[ i ] * rh[ n ]
		* A_dot[ i ][ n ];
    bias_comp_inner[ 0 ] -= temp;

    temp = 0.;
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    for( int n = 0; n < 3; ++ n )
		temp += rh[ n ] * G[ i ][ j ]
		    * W_dot[ i ][ j ][ n ];
    bias_comp_inner[ 0 ] -= 0.5 * temp;

    temp = 0.;
    for( int n = 0; n < 3; ++ n )
	temp += U_dot[ n ] * rh[ n ];
    bias_comp_inner[ 1 ] += 2. * temp;
	
    temp = 0.;
    for( int i = 0; i < 3; ++ i )
	for( int n = 0; n < 3; ++ n )
	    temp += g[ i ] * A_dot_10[ i ][ n ]
		* rh[ n ];
    bias_comp_inner[ 1 ] -= 2. * temp;

    temp = 0.;
    for( int i = 0; i < 3; ++ i )
	for( int j = 0; j < 3; ++ j )
	    for( int n = 0; n < 3; ++ n )
		temp += rh[ n ] * G[ i ][ j ]
		    * U[ i ] * A_dot[ j ][ n ];
    bias_comp_inner[ 1 ] -= 2. * temp;


    temp = 0.;
    for( int n = 0; n < 3; ++ n )
	temp += U_20_dot[ n ] * rh[ n ];
    bias_comp_inner[ 2 ] += temp;

    temp = 0.;
    for( int i = 0; i < 3; ++ i )
	for( int n = 0; n < 3; ++ n )
	    temp += g[ i ] * U[ i ]
		* U_dot[ n ] * rh[ n ];
    bias_comp_inner[ 2 ] -= 2. * temp;
    bias_comp_inner[ 3 ] -= 2. * temp;

    temp = 0.;
    for( int i = 0; i < 3; ++ i )
	for( int n = 0; n < 3; ++ n )
	    temp += g[ i ] * rh[ n ]
		* A_dot[ i ][ n ];
    bias_comp_inner[ 3 ] -= xi_R * temp;

    temp = 0.;
    for( int n = 0; n < 3; ++ n )
	temp += U_11_dot[ n ] * rh[ n ];
    bias_comp_inner[ 3 ] += temp;


    temp = 0.;
    for( int n = 0; n < 3; ++ n )
	temp += rh[ n ] * U_dot[ n ];
    bias_comp_inner[ 4 ] += 2. * xi_R * temp;
	
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

void pair_v::M1( const double & r, const double & y,
                 const double & beta )
{
    vec3 y_vec;
    y_vec.x = y * sqrt( 1 - beta * beta );
    y_vec.y = 0;
    y_vec.z = y * beta;
    M1( r, y_vec );
    return;
}

int pair_v::delta_k( const int & i, const int & j )
{
    return ( i == j ? 1 : 0 );
}

////////////////////////////////////////////////////////////
// Correlation function v12

void pair_v::v12( const double & r )
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
	    M1( r, y, beta );
	    for( int k = 0; k < num_bias_comp; ++ k )
		intg[ k ].gl_read( j, bias_comp_inner[k] );
	}
	for( int k = 0; k < num_bias_comp; ++ k )
	    intg[ k ].read( y, pow( y, 2 )
		             * intg[ k ].gl_result(  ) );

	y += dy;
    }

    for( int k = 0; k < num_bias_comp; ++ k )
	bias_comp_outer[ k ]
	    = 2 * pi * intg[ k ].result(  );

    return;
}

double pair_v::v12_L( const double & r )
{
    integral intg;
    intg.clear(  );
    k_func * kf = ( k_func * ) ( & qf->kfunc(  ));
	
    const std::vector<double> & kv = kf->kvec(  );

    double k( 0. ), jx( 0. );
    double temp_kernel( 0. );
    for( unsigned i = 0; i < kv.size(  ); ++ i )
    {
	k = kv[ i ];
	jx = k * r;
	temp_kernel = k * kf->PL_val( k )
	    * ( sin( jx ) - jx * cos( jx ) )
	    / pow( jx, 2 );
	intg.read( k, temp_kernel );
    }
    return -1. / pow( pi, 2 ) * intg.result(  );
}

void pair_v::get_v12(  )
{
    double r( 0. );
    std::cout << "Generate v12: ";
    std::cout.flush(  );

    std::vector<double> * p[ num_bias_comp ]
	={& b10b20, & b11b20, & b10b21,
	  & b12b20, & b11b21, & b10b22};

    pg.init( r_bin_num );
    const double dr = ( r_max - r_min )
	/ ( r_bin_num - 1. );
    for( int i = 0; i < r_bin_num; ++ i )
    {
	pg.show( i );
	r = dr * i + r_min;
	v12( r );
	r_buf.push_back( r );
	for( int j = 0; j < num_bias_comp; ++ j )
	    p[ j ]->push_back( bias_comp_outer[ j ] );
	v12_L_buf.push_back( v12_L( r ) );
    }

    std::cout << "Pairwise infall velocity obtained. "
	      << std::endl;
    return;
}

////////////////////////////////////////////////////////////
// Output

void pair_v::output(  )
{
    std::ofstream fout( v12_file_name.c_str(  ) );
    std::vector<double> * p[ num_bias_comp + 2 ]
	={ & r_buf, & v12_L_buf,
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


