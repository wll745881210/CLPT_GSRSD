#include "save_load.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

////////////////////////////////////////////////////////////
// Constructor, destructor, and initializer

save_load::save_load( const std::string file_path )
{
    this->file_path	= file_path;
    is_header_finished	= false;
    is_header_available = false;
    return;
}

save_load::~save_load(  )
{
    return;
}

////////////////////////////////////////////////////////////
// Read from file

void read_one_line( std::ifstream & fin )
{
    std::string line;
    std::getline( fin, line );
    if( line[ 0 ] == '#' && ( ! is_header_finished ) )
    {
	parse_header( line );
	return;
    }

    if( is_header_finished )
	
    return;
}

void save_load::parse_header( const std::string & src )
{
    std::stringstream ss;
    ss.str( src );
    
}

void save_load::
