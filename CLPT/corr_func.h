#ifndef CORR_FUNC_H
#define CORR_FUNC_H

#include "q_depend_funcs.h"
#include "lu_decomp.h"
#include "integral.h"
#include "input.h"

#include <vector>
#include <string>

struct vec3
{
    double x, y, z;
};


class corr_func
{
    ////////// Con-/destructor & initializer //////////
public:
    corr_func(  );
    ~corr_func(  );
    static void set_par( input & args );
    void initialize(  );
	
    ////////// Integration kernel //////////
protected:    // Data
    static q_func * qf;
    lu_decomp lu;
protected:    // Functions
    virtual void kernel( const double & r,
	const vec3 & y, double * bias_comp ) = 0;
    void ker_wrap( const double & r, const double & y,
	const double & beta, double * bias_comp );
    int delta_k( const int & i, const int & j );
    // Kronecker delta
	
    ////////// Correlation functions //////////
protected:    // Data
    unsigned num_bias_comp;
    static std::vector<double> rvec;
    std::vector<std::vector<double> * > corr_res;
private:    // Functions
    void calc_corr( const unsigned & i );
    virtual void post_proc(  );
public:
    void get_corr(  );

    ////////// Output //////////
public:				// Function
    void output( const std::string file_path );
	
    ////////// Mathematical const //////////
protected:
    static const double pi;
    static double y_max;
    static int    y_bin;
};

#endif

