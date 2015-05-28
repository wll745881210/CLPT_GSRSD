#include <iostream>
#include <string>

#include "q_depend_funcs.h"
#include "corr_func.h"
#include "pair_xi.h"
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

	input args( par_file_name );
	args.read(  );
	
	q_func_init q_arg;
	args.find_key( "k_file", q_arg.k_file_name, "none" );
	args.find_key( "q_file", q_arg.q_file_name, "none" );
	args.find_key( "pow_spec_file",
	               q_arg.pow_spec_name, "none" );
	q_func * p_qf = q_func::get_instance(  );
	p_qf->set_par( q_arg );

	
	corr_func::set_par( args );

	pair_xi xi;
	xi.initialize(  );
	xi.get_corr(  );
	
	std::string xi_path;
	args.find_key( "xi_file", xi_path, "xi"  );
	xi.output( xi_path );

	// pair_v v12;
	// v12.set_par( v_arg );
	// v12.get_v12(  );
	// v12.output(  );

	// pair_s s12;
	// s12.set_par( s_arg );
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

