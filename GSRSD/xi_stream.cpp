#include "integral.h"
#include "xi_stream.h"

#include <iostream>
#include <fstream>
#include <cmath>

////////////////////////////////////////////////////////////
// Constructor, destructor and initializer

xi_stream::xi_stream(  )
{

}

xi_stream::~xi_stream(  )
{
	
}

void xi_stream::set_par( const xi_stream_init & arg )
{
    this->fb11b20 = arg.fb11b20;
    this->fb10b21 = arg.fb10b21;
    this->fb11b21 = arg.fb11b21;
    this->fb12b20 = arg.fb12b20;
    this->fb10b22 = arg.fb10b22;
	
    this->y_spanning = arg.y_spanning;
    this->dy = arg.dy;
    this->f_v = arg.f_v;
    this->sigma_p_sim_100 = arg.sigma_p_sim_100;
    this->read_xi( arg.xi_file_name );
    this->read_v( arg.v_file_name );
    this->read_sigma( sigma_p, arg.sp_file_name );
    this->read_sigma( sigma_v, arg.sv_file_name );

    set_s_mu_buf( arg.s_min, arg.s_max, arg.s_bin_num,
		  arg.wedge_bin_num );
    this->out_file_name = arg.out_file_name;
    return;;
}

void xi_stream::clear(  )
{
    xi.clear(  );
    v.clear(  );
    sigma_p.clear(  );
    sigma_v.clear(  );
    s_buf.clear(  );
    xi_s_0_buf.clear(  );
    xi_s_2_buf.clear(  );
    xi_s_4_buf.clear(  );
    return;
}

////////////////////////////////////////////////////////////
// Real space functions

void xi_stream::read_to_buf
( const std::string & file_name, int col_num )
{
    std::ifstream fin( file_name.c_str(  ) );
    if( !fin )
	throw file_name.c_str(  );
    read_buf.clear(  );

    double current_val( 0. );
    while( !fin.eof(  ) )
    {
        fin >> current_val;
        read_buf.push_back( current_val );
    }
    this->col_num = col_num;
    this->row_num = read_buf.size(  ) / col_num;
    return;
}

double xi_stream::read_buf_at( int row, int col )
{
    return read_buf[ row * col_num + col ];
}

void xi_stream::read_xi( const std::string & file_name )
{
    xi.clear(  );
    xi_L.clear(  );
    read_to_buf( file_name, 8 );
    for( int i = 0; i < row_num; ++ i )
    {
        const double r = read_buf_at( i, 0 );
        const double val = read_buf_at( i, 2 )
            + fb11b20 * read_buf_at( i, 3 )
            + fb10b21 * read_buf_at( i, 4 )
            + fb12b20 * read_buf_at( i, 5 )
            + fb11b21 * read_buf_at( i, 6 )
            + fb10b22 * read_buf_at( i, 7 );
        xi.read( r, val );
        xi_L.read( r, read_buf_at( i, 1 ) );
    }	
    std::cout << file_name << " is read."<< std::endl;
    return;
}

void xi_stream::read_v( const std::string & file_name )
{
    v.clear(  );
    v_L.clear(  );
    read_to_buf( file_name, 8 );
    for( int i = 0; i < row_num; ++ i )
    {
        const double r = read_buf_at( i, 0 );
        const double val = read_buf_at( i, 2 )
            + fb11b20 * read_buf_at( i, 3 )
            + fb10b21 * read_buf_at( i, 4 )
            + fb12b20 * read_buf_at( i, 5 )
            + fb11b21 * read_buf_at( i, 6 )
            + fb10b22 * read_buf_at( i, 7 );
        v.read( r, val * f_v / ( 1. + xi( r ) ) );
        v_L.read( r, read_buf_at( i, 1 ) * f_v
                  * ( 1. + fb11b20 ) );
    }
    std::cout << file_name << " is read."<< std::endl;
    return;
}

void xi_stream::read_sigma
( interp1 & interp, const std::string & file_name )
{
    interp.clear(  );
    read_to_buf( file_name, 5 );
    for( int i = 0; i < row_num; ++ i )
    {
        const double r = read_buf_at( i, 0 );
        const double val = read_buf_at( i, 1 )
            + fb11b20 * read_buf_at( i, 2 )
            + fb10b21 * read_buf_at( i, 3 )
            + fb12b20 * read_buf_at( i, 4 );
        interp.read( r, val * f_v*f_v / ( 1. + xi( r ) ) );
    }	
    std::cout << file_name << " is read."<< std::endl;
    return;
}

////////////////////////////////////////////////////////////
// Integrand and integration

void xi_stream::set_s_mu_buf
( double s_min, double s_max, int s_bin_num, int wedge_bin_num )
{
    const double ds = ( s_max - s_min )
	/ double ( s_bin_num - 1 );
    for( int i = 0; i < s_bin_num; ++ i )
	s_buf.push_back( i * ds + s_min );
    const double dmu_wedge = 1. / double( wedge_bin_num );
    for( int i = 0; i < wedge_bin_num + 1; ++ i )
	wedge_mu_buf.push_back( i * dmu_wedge );	
    return;
}

double xi_stream::inner_integrand( const double & r_sigma,
				   const double & r_pi,
				   const double & y )
{
    const double r = sqrt( pow( r_sigma, 2 ) + pow( y, 2 ) );
    const double mu = y / r;
    if( r < 3 )
	return 0.;
	
    const double sigma2 = sigma_p( r ) * pow( mu, 2 )
	+ sigma_v( r ) * ( 1. - pow( mu, 2 ) ) * 0.5
        - sigma_shift;
    if( sigma2 < 0. )
	return 0.;
    const double exp_idx = -0.5 / sigma2
	* pow( r_pi - y - mu * v( r ), 2 );
    const double res = ( 1. + xi( r ) ) * exp( exp_idx )
	/ sqrt( 2 * pi * sigma2 );
    return res;
}

double xi_stream::inner_integration( const double & s,
				     const double & mu_s )
{
    const double r_sigma = s * sqrt( 1 - pow( mu_s, 2 ) );
    const double r_pi    = s * mu_s;
    intg.clear(  );
    double y = r_pi - y_spanning;
    double temp_integrand( 0. );
    while( y < r_pi + y_spanning )
    {
	temp_integrand = inner_integrand( r_sigma, r_pi, y );
	intg.read( y, temp_integrand );
	y += dy;
    }
    return intg.result(  ) - 1.;
}

double xi_stream::outer_integration( int order,
				     const double & s )
{
    intg.gl_clear(  );
    for( int i = 0; i < intg.gl_num; ++ i )
    {
	const double mu_s = intg.gl_xi( i );
	const double temp_xi = inner_integration( s, mu_s );
	intg.gl_read( temp_xi* legendre( order, mu_s ) );
    }
    return intg.gl_result(  ) * ( 2. * order + 1. ) * 0.5;
}

double xi_stream::wedge_integration
( const double & mu_min, const double & mu_max,
  const double & s )
{
    const double m_mu = ( mu_max - mu_min ) * 0.5;
    const double p_mu = ( mu_max + mu_min ) * 0.5;

    intg.gl_clear(  );
    for( int i = 0; i < intg.gl_num; ++ i )
    {
	const double mu_s = p_mu + intg.gl_xi( i ) * m_mu;
	const double temp_integrand
	    = inner_integration( s, mu_s );
	intg.gl_read( temp_integrand );		
    }
    return intg.gl_result(  ) * 0.5;
}

void xi_stream::gen_xi_s(  )
{
    std::cout << "Generating xi_s: ";
    std::cout.flush(  );
	
    pg.init( s_buf.size(  ) );
    xi_s_0_buf.clear(  );
    xi_s_2_buf.clear(  );
    xi_s_4_buf.clear(  );
    this->sigma_shift = sigma_p( 100. )- sigma_p_sim_100;
	
    std::vector<double> temp_wedge( wedge_mu_buf.size(  ) - 1 );
	
    for( unsigned i = 0; i < s_buf.size(  ); ++ i )
    {
	pg.show( i );
	const double & s = s_buf[ i ];
	xi_s_0_buf.push_back( outer_integration( 0, s ) );
	xi_s_2_buf.push_back( outer_integration( 2, s ) );
	xi_s_4_buf.push_back( outer_integration( 4, s ) );
		
	for( unsigned j = 0; j < wedge_mu_buf.size(  ) - 1; ++ j )
	{
	    const double mu_min = wedge_mu_buf[ j ];
	    const double mu_max = wedge_mu_buf[ j + 1 ];
	    temp_wedge[ j ]
		= wedge_integration( mu_min, mu_max, s );
	}
	xi_wedge_buf.push_back( temp_wedge );
    }
    return;
}

////////////////////////////////////////////////////////////
// Output

void xi_stream::output(  )
{
    std::ofstream fout( out_file_name.c_str(  ) );
    std::string file_name_append = out_file_name + "_appendix";
    std::ofstream fout_app( file_name_append.c_str(  ) );
    std::string file_name_wedge = out_file_name + "_wedge";
    std::ofstream fout_wedge( file_name_wedge.c_str(  ) );
    std::string file_name_2d = out_file_name + "_2d";
    std::ofstream fout_2d( file_name_2d.c_str(  ) );

    for( unsigned i = 0; i < s_buf.size(  ); ++ i )
    {
        const double s = s_buf[ i ];
        fout << s << ' ' << xi_s_0_buf[ i ] << ' '
             << xi_s_2_buf[ i ] << ' ' << xi_s_4_buf[ i ]
             << '\n';
        fout_app << s << ' ' << xi_L( s ) << ' ' << xi( s )
                 << ' ' << v_L( s ) << ' ' << v( s ) << ' '
                 << sigma_p( s ) << ' ' << sigma_v( s ) << '\n';

	fout_wedge << s << ' ';
	const std::vector<double> & wedge_row = xi_wedge_buf[ i ];
	for( unsigned j = 0; j < wedge_row.size(  ); ++ j )
	    fout_wedge << wedge_row[ j ] << ' ';
	fout_wedge << '\n';

	for( int j = intg.gl_num / 2; j < intg.gl_num; ++ j )
	{
	    const double mu = intg.gl_xi( j );
	    fout_2d << s << ' ' << mu << ' '
		    << inner_integration( s, mu ) << '\n';
	}
    }

    std::cout << "Output saved to file: \n"
	      << out_file_name << '\n'
              << file_name_append << '\n'
	      << file_name_wedge << '\n'
	      << file_name_2d << std::endl;
    return;
}

////////////////////////////////////////////////////////////
// Misc functions

double xi_stream::legendre( const int & order,
                            const double & x )
{
    switch( order )
    {
    case 0:
	return 1.;
    case 2:
	return 0.5 * ( 3. * pow( x, 2 ) - 1. );
    case 4:
	return ( 35. * pow( x, 4 ) - 30. * pow( x, 2 )
		 + 3. ) / 8.;
    default:
	throw "Inconsistent order of Legendre polynomials.";
    }
    return 0.;
}

