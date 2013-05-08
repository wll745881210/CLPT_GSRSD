#ifndef INPUT_H
#define INPUT_H

#include <fstream>
#include <string>
#include <vector>
#include "corr_func.h"
#include "q_depend_funcs.h"

class input
{
public:
	input(  )
		: fin( "par.txt" ), length( 0 )  {  };
	input( const std::string & FILE_NAME )
		: fin( FILE_NAME.c_str(  ) ), length( 0 ) {  };
	~input(  ){  };

	void read(  );
	void get_init( q_func_init & q_arg,
				   corr_func_init & c_arg,
				   corr_func_init & v_arg,
				   corr_func_init & s_arg  );

private:
	std::ifstream fin;
	std::vector<std::string> item_name;
	std::vector<std::string> value;
	int	length;
	
	void sort_item(  );
	void get_general(  );
};

#endif

