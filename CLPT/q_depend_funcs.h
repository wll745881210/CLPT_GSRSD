#ifndef Q_DEPEND_FUNCS_H
#define Q_DEPEND_FUNCS_H

#include "k_depend_funcs.h"
#include "q_depend_funcs_single.h"
#include "integral.h"
#include "prog_bar.h"
#include <vector>
#include <string>
#include <map>

////////////////////////////////////////////////////////////
// Interface for q-dependent functions 

struct q_func_init
{
    std::string pow_spec_name;
    std::string k_file_name;
    std::string q_file_name;
};

struct q_func_vals
{
    double xi_L;
    double V_112_1, V_112_3, T_112;
    double U_1,     U_3;
    double U_2_20,  U_2_11;
    double X_11,    X_22,    X_13;
    double Y_11,    Y_22,    Y_13;
    double X_12_10, Y_12_10;
};

class q_func
{
    ////////// Datatype def and declaration //////////
private:
    typedef std::vector<double> dvec;
    typedef std::map<double, double>::iterator itr_PL;
	
    ////////// Con-/destructor & initializer //////////
private:
    q_func(  );
    ~q_func(  );
    static q_func * singleton;
public:
    static q_func * get_instance(  );
    static void     del_instance(  );
    void set_par( const q_func_init & arg );

    ////////// Do all preparations //////////
private:
    void cal_all ( std::string pow_spec_name );
    void load_k  ( std::string pow_spec_name,
	           std::string k_file_name );
    void load_all( std::string pow_spec_name,
	           std::string k_file_name,
	           std::string q_file_name );

    ////////// Spherical Bessel functions //////////
private:
    double sph_bessel_j( int n, const double & x );
    // n for order.

    ////////// Various functions //////////
private:    // Data
    std::vector<q_func_single *> q_func_vec;
private:    // Function
    void set_func(  );
    void get_func_val(  );

    ////////// Function value output //////////
public:
    void var_func( const double & q, q_func_vals & res );

    ////////// Save and load //////////
private:			// Function
    void save_q_func( std::string file_name );
    void load_q_func( std::string file_name );

    ////////// Mathematical constants/func //////////
private:			// Data
    static const double nearly_0, nearly_inf;
    static const int    k_intg_points_multip;
};

#endif
