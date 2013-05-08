#include <iostream>
#include <string>

#include "q_depend_funcs.h"
#include "corr_func.h"
#include "pair_v.h"
#include "pair_s.h"
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
		q_func_init q_arg;
		corr_func_init c_arg, v_arg, s_arg;
		in.read(  );
		in.get_init( q_arg, c_arg, v_arg, s_arg );

		q_func qf;
		qf.set_par( q_arg );

		pair_s s12;
		s12.set_par( s_arg, qf );
		s12.get_s12(  );
		s12.output(  );
		
		corr_func xi;
		xi.set_par( c_arg, qf );
		xi.get_xi(  );
		xi.output(  );

		pair_v v12;
		v12.set_par( v_arg, qf );
		v12.get_v12(  );
		v12.output(  );
	}
	catch( const char * err )
	{
		std::cerr << "Error: " << err
				  << std::endl;
		return -1;
	}
	

	return 0;
}

