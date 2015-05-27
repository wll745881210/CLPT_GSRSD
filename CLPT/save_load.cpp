#include "save_load.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

////////////////////////////////////////////////////////////
// Constructor, destructor, and initializer

save_load::save_load( const std::string file_path )
{
    this->file_path	  = file_path;
    this->is_using_header = false;
    this->is_vec_assigned = false;
    return;
}

save_load::~save_load(  )
{
    for( unsigned i = 0; i < buf_in.size(  ); ++ i )
	delete ( buf_in[ i ] );
    return;
}

////////////////////////////////////////////////////////////
// Read from file

void save_load::read_one_line( std::ifstream & fin )
{
    std::string line;
    std::getline( fin, line );

    if( line[ 0 ] == '#' )
    {
	if( is_using_header )
	    header = line;
    }
    else if( ! is_vec_assigned )
	parse_head_line( line );
    else
	parse_data( line );

    return;
}

void save_load::parse_head_line( const std::string & src )
{
    std::stringstream ss;
    ss.str( src );
    if( is_using_header )
    {
	std::string name;
	while( ! ss.eof(  ) )
	{
	    ss >> name;
	    std::vector<double> * p_i
		= new std::vector<double>;

	    std::pair<std::string, std::vector<double> * >
		pair_i( name, p_i );
	    buf_in.push_back( p_i );
	    map_in.insert( pair_i );
	}
    }
    else
    {
	double val( 0. );
	while( ! ss.eof(  ) )
	{
	    ss >> val;
	    std::vector<double> * p_i
		= new std::vector<double>;
	    buf_in.push_back( p_i );
	    p_i->push_back( val );
	}
    }

    is_vec_assigned = true;
    return;
}

void save_load::parse_data( const std::string & src )
{
    std::stringstream ss;
    ss.str( src );
    double val( 0. );
    for( unsigned i = 0; i < buf_in.size(  ); ++ i )
    {
	ss >> val;
	buf_in[ i ]->push_back( val );
	if( ss.eof(  ) )
	    throw "Incorrect length of data row.";
    }
    return;
}

void save_load::set_using_header( const bool use_header )
{
    this->is_using_header = use_header;
    return;
}

void save_load::read_file(  )
{
    std::ifstream fin( file_path.c_str(  ) );
    while( ! fin.eof(  ) )
	read_one_line( fin );
    return;
}


const std::vector<double> & save_load::operator[]
( const int & i )
{
    return ( * buf_in[ i ] );
}

const std::vector<double> & save_load::operator[]
( const std::string name )
{
    std::map<std::string, std::vector<double> * >::iterator
	p = map_in.find( name );
    
    if( p == map_in.end(  ) )
	throw "Unable to find key " + name;

    return ( * p->second );
}

////////////////////////////////////////////////////////////
// Write to file

void save_load::write_one_line
( const unsigned & line_num, std::ofstream & fout )
{
    try
    {
	for( unsigned i = 0; i < buf_out.size(  ); ++ i )
	    fout << buf_out[ i ]->at( line_num ) << ' ';
	fout << '\n';
    }
    catch( ... )
    {
	std::string err = std::to_string( line_num );
	err = "Unable to write line " + err;
	throw err.c_str(  );
    }

    return;
}

void save_load::add_vec( const std::vector<double> & src,
                         const std::string name )
{
    buf_out.push_back( & src );
    if( name != "none" )
	header = header + " " + name;
    return;
}

void save_load::write(  )
{
    if( buf_out.empty(  ) )
	return;

    const unsigned l = buf_out[ 0 ]->size(  );
    for( unsigned i = 0; i < buf_out.size(  ); ++ i )
	if( buf_out[ i ]->size(  ) != l )
	    throw "Uneven length of vectors";

    std::ofstream fout( file_path.c_str(  ) );
    for( unsigned j = 0; j < l; ++ j )
	this->write_one_line( j, fout );
    return;
}
