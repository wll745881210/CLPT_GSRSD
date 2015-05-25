#include "q_depend_funcs.h"
#include "prog_bar.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>

////////////////////////////////////////////////////////////
// Static variables

const double	q_func::nearly_0	     = 1.e-3;
const double	q_func::nearly_inf	     = 2.e2;
const double	q_func::pi		     = 3.141592654;
const double	q_func::one_over_pi2	     = 0.050660592;
const int	q_func::k_intg_points_multip = 3;


////////////////////////////////////////////////////////////
// Constructor, desctructor and initializer

q_func::q_func(  )
{

}

q_func::~q_func(  )
{

}

void q_func::set_par( const q_func_init & arg )
{
    if( arg.k_file_name == "none" )
	cal_all( arg.pow_spec_name );
    else if( arg.q_file_name == "none" )
	load_k( arg.pow_spec_name, arg.k_file_name );
    else
	load_all( arg.pow_spec_name,
	    arg.k_file_name,
	    arg.q_file_name );
    return;
}


////////////////////////////////////////////////////////////
// To cal(culate) or not to cal, it is a problem.

void q_func::cal_all( std::string pow_spec_name )
{
    kf.load_PL( pow_spec_name );
    kf.get_Q_func(  );
    kf.get_R_func(  );
    kf.save_k_func( "../data/k_func.txt" );
    get_var_func(  );
    save_q_func( "../data/q_func.txt" );
    return;
}

void q_func::load_k( std::string pow_spec_name,
    std::string k_file_name )
{
    kf.load_PL( pow_spec_name );
    kf.load_k_func( k_file_name );
    get_var_func(  );
    save_q_func( "../data/q_func.txt" );
    return;
}

void q_func::load_all( std::string pow_spec_name,
    std::string k_file_name,
    std::string q_file_name )
{
    kf.load_PL( pow_spec_name );
    kf.load_k_func( k_file_name );
    load_q_func( q_file_name );
    return;
}

////////////////////////////////////////////////////////////
// Functions of k

const k_func & q_func::kfunc(  )
{
    return kf;
}

////////////////////////////////////////////////////////////
// Spherical Bessel functions

double q_func::sph_bessel_j( int n, const double & x )
{
    if( fabs( x ) < nearly_0 )
	switch( n )
	{
	case 0:
	    return 1. - x*x / 6.;
	case 1:
	    return x / 3. - pow( x, 3 ) / 30.;
	case 2:
	    return x*x / 15.;
	case 3:
	    return pow( x, 3 ) / 105.;
	default:
	    throw "Bessel function error.";
	}
    switch( n )
    {
    case 0:
	return sin( x ) / x;
    case 1:
	return sin( x ) / pow( x, 2 )
	    - cos( x ) / x;
    case 2:
	return ( 3. / pow( x, 2 ) - 1 ) * sin( x ) / x
	    - 3 * cos( x ) / pow( x, 2 );
    case 3:
	return ( 15. / pow( x, 3 ) - 6. / x ) * sin( x ) / x
	    - ( 15. / pow( x, 2 ) - 1. ) * cos( x ) / x;
    default:
	throw "Bessel function error.";
    }
    return 0.;
}

////////////////////////////////////////////////////////////
// Various functions

const std::vector<double> & q_func::qvec(  )
{
    return q_buf;
}

void q_func::get_var_func(  )
{
    double q( 0. ), k( 0. );
    std::pair<double, int> q_index_pair;
    double temp[ 3 ];

    const std::vector<double> & kv = kf.kvec(  );

    const double k_max = kv[ kv.size(  ) - 1 ]/2.;
    const double k_min = kv[ 0 ];
    const double q_max = 1. / kv[ 0 ];
    const double q_min = 1. / kv[ kv.size(  ) - 1 ];

    for( unsigned i = 0; i < kv.size(  ); ++ i )
    {
	q = log( q_min )
	    + ( log( q_max ) - log( q_min ) ) * i
	    / double( kv.size(  ) - 1 );
	q_buf.push_back( exp( q ) );
    }

    const unsigned n_k_intg
	= k_intg_points_multip * kv.size(  );
    const double d_logk
	= ( log( k_max ) - log( k_min ) )
	/ double( n_k_intg - 1 );
    for( unsigned i = 0; i < n_k_intg; ++ i )
    {
	k = exp( log( k_min ) + i * d_logk );
	k_intg_buf.push_back( k );
    }
    std::cout << "Generating q-dependent functions: ";
    pg.init( q_buf.size(  ) );

    for( unsigned i = 0; i < q_buf.size(  ); ++ i )
    {
	pg.show( i );
	q = q_buf[ i ];
	q_index_pair.first = q;
	q_index_pair.second = i;
	q_index_buf.insert( q_index_pair );

	xi_L_buf.push_back( xi_L_intg( q ) );

	temp[ 0 ] = S_112_intg( q );
	temp[ 1 ] = V_112_1_intg( q ) + temp[ 0 ];
	temp[ 2 ] = V_112_3_intg( q ) + temp[ 0 ];
	V_112_1_buf.push_back( temp[ 1 ] );
	V_112_3_buf.push_back( temp[ 2 ] );
	S_112_buf.push_back( temp[ 0 ] );
	T_112_buf.push_back( T_112_intg( q ) );

	U_1_buf.push_back( U_1_intg( q ) );
	U_3_buf.push_back( U_3_intg( q ) );
	U_2_20_buf.push_back( U_2_20_intg( q ) );
	U_2_11_buf.push_back( U_2_11_intg( q ) );

	X_11_buf.push_back( X_11_intg( q ) );
	X_22_buf.push_back( X_22_intg( q ) );
	X_13_buf.push_back( X_13_intg( q ) );
	Y_11_buf.push_back( Y_11_intg( q ) );
	Y_22_buf.push_back( Y_22_intg( q ) );
	Y_13_buf.push_back( Y_13_intg( q ) );

	X_12_10_buf.push_back( X_12_10_intg( q ) );
	Y_12_10_buf.push_back( Y_12_10_intg( q ) );
    }

    std::cout << "Various f(q) are obtained." << std::endl;
    return;
}

double q_func::xi_L_intg( const double & q )
{
    intg.clear(  );
    double k( 0. ), PL( 0. );
    double jx( 0. );	// Argument of sperical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	PL = kf.PL_val( k );
	jx = k * q;
	kernel_val = pow( k, 2 ) * PL
	    * sph_bessel_j( 0, jx );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}

double q_func::V_112_1_intg( const double & q )
{
    intg.clear(  );
    double k( 0. );
    double jx( 0. );	// Argument of sperical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	jx = k * q;
	kernel_val = -3. / 7. * kf.R_val( 1, k )
	    * sph_bessel_j( 1, jx ) / k;
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}

double q_func::V_112_3_intg( const double & q )
{
    intg.clear(  );
    double k( 0. );
    double jx( 0. );	// Argument of sperical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	jx = k * q;
	kernel_val = -3. / 7. * kf.Q_val( 1, k )
	    * sph_bessel_j( 1, jx ) / k;
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}

double q_func::S_112_intg( const double & q )
{
    intg.clear(  );
    double k( 0. );
    double jx( 0. );	// Argument of sperical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	jx = k * q;
	kernel_val = 3. / 7. / k
	    * ( 2. * kf.R_val( 1, k )
		+ 4. * kf.R_val( 2, k )
		+ kf.Q_val( 1 ,k ) + 2. * kf.Q_val( 2, k ) )
	    * sph_bessel_j( 2, jx ) / jx;
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}

double q_func::T_112_intg( const double & q )
{
    intg.clear(  );
    double k( 0. );
    double jx( 0. );	// Argument of sperical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	jx = k * q;
	kernel_val = -3. / 7. / k
	    * ( 2. * kf.R_val( 1, k )
		+ 4. * kf.R_val( 2, k )
		+ kf.Q_val( 1, k ) + 2. * kf.Q_val( 2, k) )
	    * sph_bessel_j( 3, jx );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}

double q_func::U_1_intg( const double & q )
{
    intg.clear(  );
    double k( 0. ), PL( 0. );
    double jx( 0. );	// Argument of sperical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	PL = kf.PL_val( k );
	jx = k * q;
	kernel_val = k * ( -1. ) * PL
	    * sph_bessel_j( 1, jx );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}

double q_func::U_3_intg( const double & q )
{
    intg.clear(  );
    double k( 0. );
    double jx( 0. );	// Argument of sperical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	jx = k * q;
	kernel_val = k * ( -5. / 21. ) * kf.R_val( 1, k )
	    * sph_bessel_j( 1, jx );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}

double q_func::U_2_20_intg( const double & q )
{
    intg.clear(  );
    double k( 0. );
    double jx( 0. );	// Argument of sperical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	jx = k * q;
	kernel_val = k * ( -3. / 7. )
	    * kf.Q_val( 8, k) * sph_bessel_j( 1, jx );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}

double q_func::U_2_11_intg( const double & q )
{
    intg.clear(  );
    double k( 0. );
    double jx( 0. );	// Argument of sperical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	jx = k * q;
	kernel_val = k * ( -6. / 7. )
	    * ( kf.R_val( 1, k ) + kf.R_val( 2, k ) )
	    * sph_bessel_j( 1, jx );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}

double q_func::X_11_intg( const double & q )
{
    intg.clear(  );
    double k( 0. ), PL( 0. );
    double jx( 0. );	   // Argument of sperical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	PL = kf.PL_val( k );
	jx = k * q;
	kernel_val = PL * ( 2. / 3.
	    - 2. * sph_bessel_j( 1, jx )
	    / jx );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}

double q_func::X_22_intg( const double & q )
{
    intg.clear(  );
    double k( 0. );
    double jx( 0. );	// Argument of sperical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	jx = k * q;
	kernel_val = 9. / 98. * kf.Q_val( 1, k )
	    * ( 2. / 3. - 2. * sph_bessel_j( 1, jx ) / jx );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}

double q_func::X_13_intg( const double & q )
{
    intg.clear(  );
    double k( 0. );
    double jx( 0. );	// Argument of sperical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	jx = k * q;
	kernel_val = 5. / 21. * kf.R_val( 1, k )
	    * ( 2. / 3. - 2. * sph_bessel_j( 1, jx ) / jx );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}

double q_func::Y_11_intg( const double & q )
{
    intg.clear(  );
    double k( 0. ), PL( 0. );
    double jx( 0. );	// Argument of sperical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	PL = kf.PL_val( k );
	jx = k * q;
	kernel_val = PL * ( 6. * sph_bessel_j( 1, jx ) / jx
	    - 2. * sph_bessel_j( 0, jx) );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}

double q_func::Y_22_intg( const double & q )
{
    intg.clear(  );
    double k( 0. );
    double jx( 0. );	// Argument of sperical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	jx = k * q;
	kernel_val = 9. / 98. * kf.Q_val( 1, k )
	    * ( 6. * sph_bessel_j( 1, jx ) / jx
		- 2. * sph_bessel_j( 0, jx) );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}

double q_func::Y_13_intg( const double & q )
{
    intg.clear(  );
    double k( 0. );
    double jx( 0. );	// Argument of sperical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	jx = k * q;
	kernel_val = 5. / 21. * kf.R_val( 1, k )
	    * ( 6. * sph_bessel_j( 1, jx ) / jx
		- 2. * sph_bessel_j( 0, jx) );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}

double q_func::X_12_10_intg( const double & q )
{
    intg.clear(  );
    double k( 0. );
    double jx( 0. );	// Argument of sperical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	jx = k * q;
	kernel_val = 1. / 14.
	    * ( 2. * ( kf.R_val( 1, k ) - kf.R_val( 2, k ) )
		+ 3. * kf.R_val( 1, k )
		* sph_bessel_j( 0, jx )
		- 3. * ( 3. * kf.R_val( 1, k )
		    + 4. * kf.R_val( 2, k )
		    + 2. * kf.Q_val( 5, k ) )
		* sph_bessel_j( 1, jx ) / jx );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}

double q_func::Y_12_10_intg( const double & q )
{
    intg.clear(  );
    double k( 0. );
    double jx( 0. );	// Argument of sperical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	jx = k * q;
	kernel_val = 3. / 14.
	    * ( 3. * kf.R_val( 1, k ) + 4.
		* kf.R_val( 2, k )
		+ 2. * kf.Q_val( 5, k ) )
	    * sph_bessel_j( 2, jx );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}

double q_func::interp_val( const double & q, const int & i,
    const std::vector<double> & vec )
{
    if( i < 0 || i > int( q_buf.size(  ) - 2 ) )
	return 0.;
    else if( ( q_buf[ i+1 ] - q_buf[ i ] ) / q_buf[ i ]
	< nearly_0 )
	return vec[ i ];
    else
	return vec[ i ] + ( vec[ i+1 ] - vec[ i ] )
	    / ( q_buf[ i+1 ] - q_buf[ i ] )
	    * ( q - q_buf[ i ] );
}

void q_func::var_func( const double & q, q_func_vals & res )
{
    std::map<double, int>::iterator p
	= q_index_buf.lower_bound( q );
    if( p != q_index_buf.begin(  ) )
	-- p;

    const int i = p->second;

    res.xi_L = interp_val( q, i, xi_L_buf );

    res.V_112_1 = interp_val( q, i, V_112_1_buf );
    res.V_112_3 = interp_val( q, i, V_112_3_buf );
    res.T_112 = interp_val( q, i, T_112_buf );

    res.U_1 = interp_val( q, i, U_1_buf );
    res.U_3 = interp_val( q, i, U_3_buf );
    res.U_2_20 = interp_val( q, i, U_2_20_buf );
    res.U_2_11 = interp_val( q, i, U_2_11_buf );

    res.X_11 = interp_val( q, i, X_11_buf );
    res.X_22 = interp_val( q, i, X_22_buf );
    res.X_13 = interp_val( q, i, X_13_buf );
    res.Y_11 = interp_val( q, i, Y_11_buf );
    res.Y_22 = interp_val( q, i, Y_22_buf );
    res.Y_13 = interp_val( q, i, Y_13_buf );

    res.X_12_10 = interp_val( q, i, X_12_10_buf );
    res.Y_12_10 = interp_val( q, i, Y_12_10_buf );

    return;
}

////////////////////////////////////////////////////////////
// Test output

void q_func::save_q_func( std::string file_name )
{
    std::ofstream fout( file_name.c_str(  ) );
    dvec * p[ 17 ]
	= { &q_buf, &xi_L_buf, &X_11_buf,
	    &X_22_buf, &X_13_buf, &Y_11_buf,
	    &Y_22_buf, &Y_13_buf, &U_1_buf,
	    &U_3_buf, &V_112_1_buf, &V_112_3_buf,
	    &T_112_buf, &U_2_20_buf, &U_2_11_buf,
	    &X_12_10_buf, &Y_12_10_buf };

    for( unsigned i = 0; i < q_buf.size(  ); ++ i )
    {
	for( int j = 0; j < 17; ++ j )
	    fout << p[ j ]->at( i ) << ' ';
	fout << '\n';
    }
    fout << std::endl;

    return;
}

void q_func::load_q_func( std::string file_name )
{
    std::ifstream fin( file_name.c_str(  ) );
    if( !fin )
	throw "Unable to open q function file.";

    std::pair<double, int> q_index_pair;
    dvec * p[ 17 ]
	= { &q_buf, &xi_L_buf, &X_11_buf,
	    &X_22_buf, &X_13_buf, &Y_11_buf,
	    &Y_22_buf, &Y_13_buf, &U_1_buf,
	    &U_3_buf, &V_112_1_buf, &V_112_3_buf,
	    &T_112_buf, &U_2_20_buf, &U_2_11_buf,
	    &X_12_10_buf, &Y_12_10_buf };
    for( int i = 0; i < 17; ++ i )
	p[ i ]->clear(  );

    double temp( 0. );
    while( !fin.eof(  ) )
    {

	for( int i = 0; i < 17; ++ i )
	{
	    fin >> temp;
	    p[ i ]->push_back( temp );
	}
    }
    for( int i = 0; i < 17; ++ i )
	p[ i ]->pop_back(  );
    for( unsigned i = 0; i < q_buf.size(  ); ++ i )
    {
	q_index_pair.first = q_buf[ i ];
	q_index_pair.second = i;
	q_index_buf.insert( q_index_pair );
    }

    std::cout << "q functions loaded from: "
	      << file_name << std::endl;
    return;
}
