#include <iostream>
#include <string>

#include "xi_stream.h"
#include "input.h"

void driver( const std::string & par_file_name )
{
    xi_stream_init arg;
    xi_stream test;
    
    input par( par_file_name );
    par.read(  );
    arg.s_max = 130;
    arg.s_min = 1;
    arg.s_bin_num = 50;
    arg.y_spanning = 200;
    arg.dy = 0.5;
    arg.f_v = 0.74429;
    arg.sigma_p_sim_100 = 27.;
    arg.xi_file_name  = "xi.txt";
    arg.v_file_name   = "v12.txt";
    arg.sp_file_name  = "s12.txt_p";
    arg.sv_file_name  = "s12.txt_v";
    arg.out_file_name = "xi_s.txt";
    par.find_key( "fb11b20", arg.fb11b20 );
    par.find_key( "fb10b21", arg.fb10b21 );
    par.find_key( "fb11b21", arg.fb11b21 );
    par.find_key( "fb12b20", arg.fb12b20 );
    par.find_key( "fb10b22", arg.fb10b22 );
    par.find_key( "s_max", arg.s_max );
    par.find_key( "s_min", arg.s_min );
    par.find_key( "s_bin", arg.s_bin_num );
    par.find_key( "mu_bin", arg.mu_bin_num );
    par.find_key( "y_spanning", arg.y_spanning );
    par.find_key( "dy", arg.dy );
    par.find_key( "f_v", arg.f_v );
    par.find_key( "sigma_p_sim_100", arg.sigma_p_sim_100 );
    par.find_key( "xi_file_name", arg.xi_file_name );
    par.find_key( "v_file_name" , arg.v_file_name  );
    par.find_key( "sp_file_name", arg.sp_file_name );
    par.find_key( "sv_file_name", arg.sv_file_name );
    par.find_key( "out_file_name",arg.out_file_name );
        
    test.clear(  );
    test.set_par( arg );
    test.gen_xi_s(  );
    test.output(  );

    return;
}

int main( int argn, char * argv[  ] )
{
	try
	{
        std::string par_file_name;
        if( argn == 1 )
            par_file_name = "par.txt";
        else if( argn == 2 )
            par_file_name = argv[ 1 ];
        else
            throw "Incorrect specification of input file.";

        driver( par_file_name );
	}
	catch( const char * err )
	{
		std::cerr << "\nError: " << err
				  << std::endl;
	}
	return 0;
}


