#include <iostream>
#include <string>
#include <omp.h>

#include "q_depend_funcs.h"
#include "corr_func.h"
#include "pair_xi.h"
#include "pair_v.h"
#include "pair_s.h"
#include "input.h"

void corr_driver( corr_func & corr,
		  std::string specifier, input & args )
{
    corr.initialize(  );
    corr.get_corr(  );
    std::string path;
    args.find_key( specifier + "_file", path, specifier );
    corr.output( path );
    return;
}

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

	int n_thread( 0 );
	args.find_key( "n_thread", n_thread, -1 );
	if( n_thread > 0 )
	    omp_set_num_threads( n_thread );
	
	q_func * p_qf = q_func::get_instance(  );
	p_qf->initialize( args );
	
	corr_func::set_par( args );
	pair_xi xi;
	pair_v  v12;
	pair_s  s12;
	corr_driver( xi,  "xi",  args );
	corr_driver( v12, "v12", args );
	corr_driver( s12, "s12", args );
    }
    catch( const char * err )
    {
	std::cerr << "\nError: " << err << std::endl;
	return -1;
    }


    return 0;
}
