#ifndef PAIR_V_H
#define PAIR_V_H

#include "q_depend_funcs.h"
#include "corr_func.h"
#include "lu_decomp.h"
#include "integral.h"
#include "prog_bar.h"

#include <vector>
#include <string>

class pair_v
{
    ////////// Con-/destructor & initializer //////////
public:
    pair_v(  );
    ~pair_v(  );
    void set_par( const corr_func_init & v_arg,
	const q_func & q );
	
    ////////// Integration kernel //////////
private:    // Data
    q_func * qf;
    lu_decomp lu;
private:    // Functions
    void M1( const double & r, const vec3   & y );
    void M1( const double & r, const double & y,
             const double & beta );
    int delta_k( const int & i, const int & j );
	
    ////////// Corr func v12 //////////
private:    // Data
    double r_max, r_min;
    int r_bin_num;
    static const int num_bias_comp = 6;
    double bias_comp_inner[ num_bias_comp ];
    double bias_comp_outer[ num_bias_comp ];
    std::vector<double> r_buf;
    std::vector<double> v12_L_buf;
    std::vector<double> b10b20, b11b20, b10b21;
    std::vector<double> b12b20, b11b21, b10b22;
private:    // Functions
    void   v12  ( const double & r );
    double v12_L( const double & r );
public:
    void get_v12(  );

    ////////// Output //////////
private:
    std::string v12_file_name;
public:
    void output(  );
	
    ////////// Progress bar //////////
private:
    prog_bar pg;
	
    ////////// Mathematical const //////////
private:
    static const double pi;
    static const double max_y;
    static const int    num_y;
};

#endif

