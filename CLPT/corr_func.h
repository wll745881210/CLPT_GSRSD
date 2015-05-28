#ifndef CORR_FUNC_H
#define CORR_FUNC_H

#include "q_depend_funcs.h"
#include "lu_decomp.h"
#include "integral.h"
#include "prog_bar.h"
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
private:    // Data
    static q_func * qf;
    lu_decomp lu;
private:    // Functions
    virtual void kernel
    ( const double & r, const vec3 & y ) = 0;
    void kernel( const double & r, const double & y,
	         const double & beta );
    int delta_k( const int & i, const int & j );
    // Kronecker delta
	
    ////////// Correlation functions //////////
private:    // Data
    double * bias_comp;
    static std::vector<double> rvec;
    std::vector<std::vector<double> * > corr_res;
private:    // Functions
    void calc_corr( const double & r );
public:
    void get_corr(  );

    ////////// Output //////////
public:				// Function
    void output( const std::string file_path;  );
	
    ////////// Mathematical const //////////
private:
    static const double pi;
    static double y_max;
    static int    y_bin;
};

#endif

