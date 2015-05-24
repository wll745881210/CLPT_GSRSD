#ifndef Q_DEPEND_FUNCS_H
#define Q_DEPEND_FUNCS_H

#include "k_depend_funcs.h"
#include "integral.h"
#include "prog_bar.h"
#include <vector>
#include <string>
#include <map>

// All q-dependent functions in Jordan Carlson's notes

struct q_func_init
{
    std::string pow_spec_name;
    std::string k_file_name;
    std::string q_file_name;
};

struct q_func_vals
{
    double xi_L;
    double V_112_1, V_112_3;
    double T_112;
    double U_1, U_3;
    double U_2_20, U_2_11;
    double X_11, X_22, X_13;
    double Y_11, Y_22, Y_13;
    double X_12_10, Y_12_10;
};


class q_func
{
    ////////// Datatype def and declaration //////////
private:
    typedef std::vector<double> dvec;
    typedef std::map<double, double>::iterator itr_PL;
	
    ////////// Con-/destructor & initializer //////////
public:
    q_func(  );
    ~q_func(  );
    void set_par( const q_func_init & arg );

    ////////// Do all preparations //////////
private:
    void cal_all( std::string pow_spec_name );
    void load_k( std::string pow_spec_name,
	std::string k_file_name );
    void load_all( std::string pow_spec_name,
	std::string k_file_name,
	std::string q_file_name );

    ////////// Functions of k //////////
public:
    const k_func & kfunc(  );
private:
    k_func kf;
	
    ////////// Spherical Bessel functions //////////
private:
    double sph_bessel_j( int n, const double & x );
    // n for order.

    ////////// Various functions //////////
public:
    const std::vector<double> & qvec(  );
private:						// Data
    std::map<double, int> q_index_buf;
    std::vector<double> q_buf;
    std::vector<double> k_intg_buf;
    std::vector<double> xi_L_buf;
    std::vector<double> V_112_1_buf, V_112_3_buf;
    std::vector<double> S_112_buf, T_112_buf;
    std::vector<double> U_1_buf, U_3_buf;
    std::vector<double> U_2_20_buf, U_2_11_buf;
    std::vector<double> X_11_buf, X_22_buf, X_13_buf;
    std::vector<double> Y_11_buf, Y_22_buf, Y_13_buf;
    std::vector<double> X_12_10_buf, Y_12_10_buf;
private:						// Function
    void get_var_func(  );
    double xi_L_intg( const double & q );
    double V_112_1_intg( const double & q );
    double V_112_3_intg( const double & q );
    double S_112_intg( const double & q );
    double T_112_intg( const double & q );
    double U_1_intg( const double & q );
    double U_3_intg( const double & q );
    double U_2_20_intg( const double & q );
    double U_2_11_intg( const double & q );
    double X_11_intg( const double & q );
    double X_22_intg( const double & q );
    double X_13_intg( const double & q );
    double Y_11_intg( const double & q );
    double Y_22_intg( const double & q );
    double Y_13_intg( const double & q );
    double X_12_10_intg( const double & q );
    double Y_12_10_intg( const double & q );
    // Linear interpolation 
    double interp_val( const double & q, const int & i,
	const std::vector<double> & vec );
public:
    void var_func( const double & q, q_func_vals & res );

    ////////// Save and load //////////
private:
    void save_q_func( std::string file_name );
    void load_q_func( std::string file_name );

    ////////// Progress bar //////////
private:
    prog_bar pg;
	
    ////////// Mathematical constants/func //////////
private:
    static const double nearly_0 = 1.e-3;
    static const double nearly_inf = 2.e2;
    static const double pi = 3.14159265358979323846;
    static const double one_over_pi2 = 0.05066059182116889;
    static const int k_intg_points_multip = 3;
    integral intg;
};

#endif

