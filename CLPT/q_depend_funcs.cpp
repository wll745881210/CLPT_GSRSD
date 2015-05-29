#include "q_depend_funcs_single.h"
#include "q_depend_funcs.h"
#include "input.h"
#include "save_load.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>

////////////////////////////////////////////////////////////
// Static variables

const double q_func::nearly_0            ( 1.e-3   );
const double q_func::nearly_inf          ( 2.e2    );
const int    q_func::k_intg_points_multip( 3       );
q_func *     q_func::singleton           ( nullptr );

////////////////////////////////////////////////////////////
// Constructor, desctructor and initializer

q_func::q_func(  )
{
    this->set_func(  );
    return;
}

q_func::~q_func(  )
{
    return;
}

q_func * q_func::get_instance(  )
{
    if( singleton == nullptr )
	singleton = new q_func;
    return singleton;
}

void q_func::del_instance(  )
{
    delete singleton;
    singleton = nullptr;
    return;
}

void q_func::initialize( input & args )
{
    std::string k_path_in,  q_path_in,  ps_path;
    std::string k_path_out, q_path_out;

    args.find_key( "k_input",      k_path_in,  "none" );
    args.find_key( "q_input",      q_path_in,  "none" );
    args.find_key( "pow_spec_file",   ps_path, "none" );

    const std::string base( "../data/" );
    args.find_key( "k_output", k_path_out,
	           base + "k_func.txt" );
    args.find_key( "q_output", q_path_out,
	           base + "q_func.txt" );
    
    if( k_path_in == "none" )
	calc_all( ps_path, k_path_out, q_path_out );
    else if( q_path_in == "none" )
	load_k  ( ps_path, k_path_in,  q_path_out );
    else
	load_all( ps_path, k_path_in,  q_path_in  );
    return;
}


////////////////////////////////////////////////////////////
// To cal(culate) or not to cal, it is a problem.

void q_func::calc_all( const std::string pow_spec_path,
                       const std::string k_path,
                       const std::string q_path )
{
    k_func * p_kf = k_func::get_instance(  );
    p_kf->load_PL     ( pow_spec_path );
    p_kf->get_Q_func  (               );
    p_kf->get_R_func  (               );
    p_kf->save_k_func ( k_path        );
    this->get_func_val(               );
    this->save_q_func ( q_path        );
    return;
}

void q_func::load_k( const std::string pow_spec_path,
                     const std::string k_path,
                     const std::string q_path )
{
    k_func * p_kf = k_func::get_instance(  );
    p_kf->load_PL     ( pow_spec_path );
    p_kf->load_k_func ( k_path        );
    this->get_func_val(               );
    this->save_q_func ( q_path        );
    return;
}

void q_func::load_all( const std::string pow_spec_path,
                       const std::string k_path,
                       const std::string q_path )
{
    k_func * p_kf = k_func::get_instance(  );
    p_kf->load_PL    ( pow_spec_path );
    p_kf->load_k_func( k_path        );
    this->load_q_func( q_path        );
    return;
}

////////////////////////////////////////////////////////////
// Various functions

void q_func::set_func(  )
{
    class xi_L : public q_func_single
    {
    	double kernel( const double & k, const double & jx,
    	               const k_func & kf )
    	{
    	    return pow( k, 2 ) * kf.PL_val( k )
    		 * sph_bessel_j( 0, jx );
    	}
    };
    q_func_vec.push_back( new xi_L );

    class X_11 : public q_func_single
    {
	double kernel( const double & k, const double & jx,
	               const k_func & kf )
	{
	    return kf.PL_val( k ) * ( 2. / 3.
		 - 2. * sph_bessel_j( 1, jx ) / jx );
	}
    };
    q_func_vec.push_back( new X_11 );

    class X_22 : public q_func_single
    {
	double kernel( const double & k, const double & jx,
	               const k_func & kf )
	{
	    return 9. / 98. * kf.Q_val( 1, k )
		* ( 2. / 3. - 2. * sph_bessel_j( 1, jx )
		    / jx );
	}
    };
    q_func_vec.push_back( new X_22 );

    class X_13 : public q_func_single
    {
	double kernel( const double & k, const double & jx,
	               const k_func & kf )
	{
	    return 5. / 21. * kf.R_val( 1, k )
		* ( 2. / 3. - 2. * sph_bessel_j( 1, jx )
		    / jx );
	}
    };
    q_func_vec.push_back( new X_13 );
    
    class Y_11 : public q_func_single
    {
	double kernel( const double & k, const double & jx,
	               const k_func & kf )
	{
	    return kf.PL_val( k )
		* ( 6. * sph_bessel_j( 1, jx ) / jx
		    - 2. * sph_bessel_j( 0, jx ) );
	}
    };
    q_func_vec.push_back( new Y_11 );

    class Y_22 : public q_func_single
    {
	double kernel( const double & k, const double & jx,
	               const k_func & kf )
	{
	    return  9. / 98. * kf.Q_val( 1, k )
	    * ( 6. * sph_bessel_j( 1, jx ) / jx
		- 2. * sph_bessel_j( 0, jx ) );
	}
    };
    q_func_vec.push_back( new Y_22 );

    class Y_13 : public q_func_single
    {
	double kernel( const double & k, const double & jx,
	               const k_func & kf )
	{
	    return  5. / 21. * kf.R_val( 1, k )
		* ( 6. * sph_bessel_j( 1, jx ) / jx
		    - 2. * sph_bessel_j( 0, jx ) );
	}
    };
    q_func_vec.push_back( new Y_13 );

    class U_1 : public q_func_single
    {
	double kernel( const double & k, const double & jx,
	               const k_func & kf )
	{
	    return -k * kf.PL_val( k )
		* sph_bessel_j( 1, jx );
	}
    };
    q_func_vec.push_back( new U_1 );

    class U_3 : public q_func_single
    {
	double kernel( const double & k, const double & jx,
	               const k_func & kf )
	{
	    return k * ( -5. / 21. ) * kf.R_val( 1, k )
 		 * sph_bessel_j( 1, jx );
	}
    };
    q_func_vec.push_back( new U_3 );

    class V_112_1 : public q_func_single
    {
	double kernel( const double & k, const double & jx,
	               const k_func & kf )
	{
	    return -3. / 7. * kf.R_val( 1, k )
		 * sph_bessel_j( 1, jx ) / k
		+ 3. / 7. / k
		* ( 2. * kf.R_val( 1, k )
		    + 4. * kf.R_val( 2, k )
		    + kf.Q_val( 1 ,k )
		    + 2. * kf.Q_val( 2, k ) )
		* sph_bessel_j( 2, jx ) / jx;
	}
    };
    q_func_vec.push_back( new V_112_1 );

    class V_112_3 : public q_func_single
    {
	double kernel( const double & k, const double & jx,
	               const k_func & kf )
	{
	    return  -3. / 7. * kf.Q_val( 1, k )
		 * sph_bessel_j( 1, jx ) / k
		+ 3. / 7. / k
		* ( 2. * kf.R_val( 1, k )
		    + 4. * kf.R_val( 2, k )
		    + kf.Q_val( 1 ,k )
		    + 2. * kf.Q_val( 2, k ) )
		* sph_bessel_j( 2, jx ) / jx;
	}
    };
    q_func_vec.push_back( new V_112_3 );
    
    class T_112 : public q_func_single
    {
	double kernel( const double & k, const double & jx,
	               const k_func & kf )
	{
	    return -3. / 7. / k
		* ( 2. * kf.R_val( 1, k )
		    + 4. * kf.R_val( 2, k )
		    + kf.Q_val( 1, k )
		    + 2. * kf.Q_val( 2, k ) )
		* sph_bessel_j( 3, jx );
	}
    };
    q_func_vec.push_back( new T_112 );
    
    class U_2_20 : public q_func_single
    {
	double kernel( const double & k, const double & jx,
	               const k_func & kf )
	{
	    return k * ( -3. / 7. )
		 * kf.Q_val( 8, k) * sph_bessel_j( 1, jx );
	}
    };
    q_func_vec.push_back( new U_2_20 );

    class U_2_11 : public q_func_single
    {
	double kernel( const double & k, const double & jx,
	               const k_func & kf )
	{
	    return k * ( -6. / 7. )
		 * ( kf.R_val( 1, k ) + kf.R_val( 2, k ) )
		 * sph_bessel_j( 1, jx );
	}
    };
    q_func_vec.push_back( new U_2_11 );

    class X_12_10 : public q_func_single
    {
	double kernel( const double & k, const double & jx,
	               const k_func & kf )
	{
	    return  1. / 14.
		 * ( 2. * ( kf.R_val( 1, k )
			 - kf.R_val( 2, k ) )
		     + 3. * kf.R_val( 1, k )
		     * sph_bessel_j( 0, jx )
		     - 3. * ( 3. * kf.R_val( 1, k )
			    + 4. * kf.R_val( 2, k )
			    + 2. * kf.Q_val( 5, k ) )
		     * sph_bessel_j( 1, jx ) / jx );
	}
    };
    q_func_vec.push_back( new X_12_10 );

    class Y_12_10 : public q_func_single
    {
	double kernel( const double & k, const double & jx,
	               const k_func & kf )
	{
	    return 3. / 14.
		 * ( 3. * kf.R_val( 1, k ) + 4.
		     * kf.R_val( 2, k )
		     + 2. * kf.Q_val( 5, k ) )
		* sph_bessel_j( 2, jx );
	}
    };
    q_func_vec.push_back( new Y_12_10 );

    return;
}

void q_func::get_func_val(  )
{
    k_func * p_kf = k_func::get_instance(  );
    const std::vector<double> & kv = p_kf->kvec(  );

    //////////////////////////////////////////////////
    // Assigning q_buf and k_buf. Quite empirical.
    std::vector<double> q_buf, k_buf;

    const double q_max  = 1. / kv[ 0 ];
    const double q_min  = 1. / kv[ kv.size(  ) - 1 ];
    const double d_logq = ( log( q_max ) - log( q_min ) )
	                / double( kv.size(  ) - 1 );
    for( unsigned i = 0; i < kv.size(  ); ++ i )
	q_buf.push_back( exp( log( q_min ) + i * d_logq ) );

    const double k_max  = kv[ kv.size(  ) - 1 ] / 2.;
    const double k_min  = kv[ 0 ];
    const int    n_k    = k_intg_points_multip * kv.size();
    const double d_logk	= ( log( k_max ) - log( k_min ) )
	                / double( n_k - 1 );
    for( int i = 0; i < n_k; ++ i )
	k_buf.push_back( exp( log( k_min ) + i * d_logk ) );
    //////////////////////////////////////////////////

    
    std::cout << "Generating q-dependent functions: "
	      << std::flush;
    q_func_single::set_qvec( q_buf );
    q_func_single::set_kvec( k_buf );

    for( unsigned i = 0; i < q_func_vec.size(  ); ++ i )
	q_func_vec[ i ]->integrate(  );
    
    std::cout << "Done." << std::endl;
    
    return;
}

////////////////////////////////////////////////////////////
// Function value output

void q_func::var_func( const double & q, q_func_vals & res )
{
    // Ugly, huh? Okay, I delibrately avoid using Boost...
    double * p_res[ 16 ]
	= { & res.xi_L,    & res.X_11,    & res.X_22,
	    & res.X_13,    & res.Y_11,    & res.Y_22,
	    & res.Y_13,    & res.U_1,     & res.U_3,
	    & res.V_112_1, & res.V_112_3, & res.T_112,
	    & res.U_2_20,  & res.U_2_11,  & res.X_12_10,
	    & res.Y_12_10 };

    const unsigned idx = q_func_single::eval_all_idx( q );
    for( unsigned i = 0; i < q_func_vec.size(  ); ++ i )
	( * p_res[ i ] )
	    = q_func_vec[ i ]->get_val( q, idx );

    return;
}

////////////////////////////////////////////////////////////
// Save/load

void q_func::save_q_func( const std::string & file_name )
{
    save_load s( file_name );
    s.add_vec( q_func_single::get_qvec(  ) );
    for( unsigned i = 0; i < q_func_vec.size(  ); ++ i )
    	s.add_vec( q_func_vec[ i ]->get_valvec(  ) );

    s.write(  );
    return;
}

void q_func::load_q_func( const std::string & file_name )
{
    save_load l( file_name );
    l.read_file(  );

    q_func_single::set_qvec( l[ 0 ] );
    for( unsigned i = 0; i < q_func_vec.size(  ); ++ i )
	q_func_vec[ i ]->set_valvec( l[ i + 1 ] );

    std::cout << "q functions loaded from: "
    	      << file_name << std::endl;
    return;
}
