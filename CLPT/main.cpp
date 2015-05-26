#include <iostream>
#include <string>

#include "q_depend_funcs.h"
// #include "corr_func.h"
// #include "pair_v.h"
// #include "pair_s.h"
#include "input.h"

int main( int argn, char * argv[  ] )
{
    try
    {
	std::string par_file_name;
	if( argn > 2 )
	    throw "Incorrect input parameter.\n"
		"See README for usage.";
	else if( argn < 2 )
	    par_file_name = "par.txt";
	else
	    par_file_name = argv[ 1 ];

	input in( par_file_name );
	in.read(  );
	
	q_func_init q_arg;
	in.find_key( "k_file", q_arg.k_file_name, "none" );
	in.find_key( "q_file", q_arg.q_file_name, "none" );
	in.find_key( "pow_spec_file",
                     q_arg.pow_spec_name, "none" );
	q_func * p_qf = q_func::get_instance(  );
	p_qf->set_par( q_arg );

	// corr_func_init c_arg, v_arg, s_arg;
	// in.find_key( "r_max",     c_arg.r_max,     130 );
	// in.find_key( "r_min",     c_arg.r_min,     1   );
	// in.find_key( "r_bin_num", c_arg.r_bin_num, 80  );

	// v_arg = c_arg;
	// s_arg = c_arg;
	
	// in.find_key( "xi_file",  c_arg.file_name, "xi"  );
	// in.find_key( "v12_file", v_arg.file_name, "v12" );
	// in.find_key( "s12_file", s_arg.file_name, "s12" );
		
	// corr_func xi;
	// xi.set_par( c_arg, * p_qf );
	// xi.get_xi(  );
	// xi.output(  );

	// pair_v v12;
	// v12.set_par( v_arg, * p_qf );
	// v12.get_v12(  );
	// v12.output(  );

	// pair_s s12;
	// s12.set_par( s_arg, * p_qf );
	// s12.get_s12(  );
	// s12.output(  );

    }
    catch( const char * err )
    {
	std::cerr << "\nError: " << err << std::endl;
	return -1;
    }
	

    return 0;
}

