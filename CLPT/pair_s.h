#ifndef PAIR_S_H
#define PAIR_S_H

#include "q_depend_funcs.h"
#include "corr_func.h"
#include "lu_decomp.h"
#include "integral.h"
#include "prog_bar.h"

#include <vector>
#include <string>

class pair_s
{
    ////////// Con-/destructor & initializer //////////
public:
    pair_s(  );
    ~pair_s(  );
    void set_par( const corr_func_init & s_arg );
	
    ////////// Integration kernel //////////
private:    // Data
    q_func * qf;
    lu_decomp lu;
private:    // Functions
    void M2( const double & r, const vec3 & y );
    void M2( const double & r, const double & y,
	     const double & beta );
    int delta_k( const int & i, const int & j );
	
    ////////// Corr func v12 //////////
private:    // Data
    double r_max, r_min;
    int r_bin_num;
    static const int num_bias_comp = 4;
    double bias_comp_inner_p[ num_bias_comp ];
    double bias_comp_inner_v[ num_bias_comp ];
    double bias_comp_outer_p[ num_bias_comp ];
    double bias_comp_outer_v[ num_bias_comp ];
    std::vector<double> r_buf;
    std::vector<double> b0p, b11p, b21p, b12p;
    std::vector<double> b0v, b11v, b21v, b12v;
private:    // Functions
    void s12( const double & r );
public:
    void get_s12(  );

    ////////// Output //////////
private:
    std::string s12_file_name;
public:
    void output(  );
	
    ////////// Progress bar //////////
private:
    prog_bar pg;
	
    ////////// Mathematical const //////////
private:
    static const double pi;
    static const double nearly_0;
    static const double max_y;
    static const int    num_y;
};

#endif

