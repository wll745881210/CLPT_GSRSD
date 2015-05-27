#ifndef SAVE_LOAD_H
#define SAVE_LOAD_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

class save_load
{
    ////////// Con-/Destructor and initializer //////////
public:
    save_load( const std::string file_path );
    ~save_load(  );

    ////////// General data //////////
private:			// Data
    std::string file_path;
    std::vector<std::string> header;
    std::map<std::string, std::vector<double> *> buf_map;
    std::vector<std::vector<double> *>           buf_vec;

    ////////// Read from file //////////
private:			// Data
    bool is_header_available;
    bool is_header_finished;
private:			// Function
    void read_one_line( std::ifstream & fin );
    void set_
    void parse_header ( const std::string & src );
    void parse_data   ( const std::string & src );
public:				// Function
    void parse_header(  )
    void read_file(  );

    ////////// 
    
    
};

#endif
