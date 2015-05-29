#ifndef SAVE_LOAD_H
#define SAVE_LOAD_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>

class save_load
{
    ////////// Con-/Destructor and initializer //////////
public:
    save_load( const std::string file_path );
    ~save_load(  );

    ////////// General data //////////
private:			// Data
    std::string file_path;
    std::string header;

    ////////// Read from file //////////
private:			// Data
    bool is_using_header;
    bool is_vec_assigned;
    std::unordered_map <std::string, std::vector<double> *>
                       map_in;
    std::vector<std::vector<double> *>           buf_in;
private:			// Function
    void read_one_line( std::ifstream & fin );
    void parse_head_line( const std::string & src );
    void parse_data     ( const std::string & src );
public:				// Function
    void set_using_header( const bool use_header );
    void read_file(  );
    const std::vector<double> & operator[]
    ( const int & i );
    const std::vector<double> & operator[]
    ( const std::string name );

    ////////// Write to file //////////
private:			// Data
    std::vector<const std::vector<double> * > buf_out;
private:			// Function
    
    void write_one_line( const unsigned & line_num,
	                 std::ofstream & fout );
public:				// Function
    void add_vec( const std::vector<double> & src,
	          const std::string name = "none" );
    void write(  );
};

#endif
