#include <iostream>
#include <string>

#include "xi_stream.h"
#include "input.h"

void driver( const std::string & par_file_name )
{
    xi_stream xi_s;
    
    input arg( par_file_name );
    arg.read(  );
    xi_s.clear(  );
    xi_s.set_par( arg );
    xi_s.gen_xi_s(  );
    xi_s.output(  );

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
	std::cerr << "\nError: " << err << std::endl;
    }
    return 0;
}


