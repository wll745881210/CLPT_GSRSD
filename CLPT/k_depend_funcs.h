#ifndef K_DEPEND_FUNCS_H
#define K_DEPEND_FUNCS_H

#include <map>
#include <vector>
#include <string>
#include "integral.h"
#include "prog_bar.h"

class k_func
{
    ////////// Type definition //////////
private:
    typedef std::map<double, int>::iterator itr_idx_map;
    typedef std::vector<double> dvec;
	
    ////////// Con-/destructor //////////
private:
    k_func(  );
    ~k_func(  );
    static k_func * singleton;
public:
    static k_func * get_instance(  );
    static void     del_instance(  );
	
    ////////// Interpolation //////////
private:    // Data
    std::map<double, int> idx_map;
private:    // Function
    double interp( const double & k,
	           const std::vector<double> & vec );

    ////////// Linear power spectrum //////////
public:
    void load_PL( std::string file_name );
    double PL_val( const double & k );
    const std::vector<double> & kvec(  );
private:    // Data
    double k_min, k_max;
    std::vector<double> k_buf;
    std::vector<double> PL_buf;
	
    ////////// Q functions //////////
public:
    double Q_val( const int n, const double & k );
    void get_Q_func(  );
private:    // Data
    std::vector<double> Q_buf_1, Q_buf_2, Q_buf_3;
    std::vector<double> Q_buf_4, Q_buf_5, Q_buf_6;
    std::vector<double> Q_buf_7, Q_buf_8;
private:    // Function
    double Q_kernel( int n, const double & r,
	const double & x );
    double Q_inner_integration( int n, const double & r,
	const double & k );
    double Q_outer_integration( int n, const double & k );

    ////////// R functions //////////
public:
    double R_val( const int n, const double & k );
    void get_R_func(  );
private:    // Data
    std::vector<double> R_buf_1, R_buf_2;
private:    // Function
    double R_inner_integration( int n, const double & r );
    double R_outer_integration( int n, const double & k );
	
    ////////// Save and load //////////
public:
    void save_k_func( std::string file_name );
    void load_k_func( std::string file_name );
		
    ////////// Progress bar //////////
private:
    prog_bar pg;

    ////////// Mathematical constants/func //////////
private:
    static const double nearly_0, nearly_inf;
    static const double pi,       one_over_pi2;
};



#endif
