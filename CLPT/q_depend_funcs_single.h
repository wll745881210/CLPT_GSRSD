#ifndef Q_DEPEND_FUNCS_SINGLE_H
#define Q_DEPEND_FUNCS_SINGLE_H

#include "k_depend_funcs.h"
#include <vector>
#include <string>
#include <map>

////////////////////////////////////////////////////////////
// All q-dependent functions in Jordan Carlson's notes

////////////////////////////////////////////////////////////
// Virtual base class for individual function.
// I avoid lambdas so that it can be compiled on those
// systems without C++11 support.

class q_func_single
{
    ////////// Constructor/Destructor //////////
public:				// Function
     q_func_single(  );
    ~q_func_single(  );

    ////////// Integration kernel specification //////////
private:    			// Function
    virtual double kernel
    ( const double & k, const double & jx,
      const k_func & kf ) = 0;

    ////////// Integrator //////////
private:			// Data
    k_func * p_kf;
protected:			// Function
    static double sph_bessel_j
    ( const int & n, const double & x );
    double integrate_single( const double & q );
public:				// Function
    void integrate(  );

    ////////// Integration results //////////
private:			// Data
    static std::vector<double>   qvec, kvec;
    static std::map<double, int> q_idx_buf;
    std::vector<double> valvec;
public:				// Function
    static void set_kvec
    ( const std::vector<double> & k_src );
    static void set_qvec
    ( const std::vector<double> & q_src );
    static const std::vector<double> & get_qvec(  );
    double & valvec_at( const int & i );

    ////////// Interpolation //////////
private:			// Data
    static int idx_current;
    static double q_current;
public:				// Function-like
    static void eval_all_idx( const double & q );
    double get_val(  );
};

#endif

