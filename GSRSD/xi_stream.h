#ifndef XI_STREAM_H
#define XI_STREAM_H

#include "integral.h"
#include "prog_bar.h"
#include "interp1.h"

#include <string>
#include <vector>
#include <map>

struct xi_stream_init
{
	double s_max, s_min;
	int s_bin_num, mu_bin_num;
	double y_spanning, dy;
    double fb11b20, fb10b21, fb11b21;
    double fb12b20, fb10b22;
    double f_v;
    double sigma_p_sim_100;
	std::string xi_file_name;
	std::string v_file_name;
	std::string sp_file_name;
	std::string sv_file_name;	
	std::string out_file_name;
};

class xi_stream
{
	////////// Con-/destructor & initializer //////////
public:
	xi_stream(  );
	~xi_stream(  );
	void set_par( const xi_stream_init & arg );
	void clear(  );

	////////// Real space func //////////
private:						// Data
	interp1 xi, xi_L, v, v_L, sigma_p, sigma_v;
    double fb11b20, fb10b21, fb11b21;
    double fb12b20, fb10b22;
    double f_v;
    double sigma_p_sim_100, sigma_shift;
    std::vector<double> read_buf;
    int col_num, row_num;
private:						// Function
    void read_to_buf( const std::string & file_name,
                      int col_num );
    inline double read_buf_at( int col, int row );
	void read_xi( const std::string & file_name );
	void read_v(  const std::string & file_name );
	void read_sigma( interp1 & interp,
                     const std::string & file_name );

	////////// Integration/integrand //////////
private:						// Data
	double y_spanning, dy;
	std::vector<double> s_buf, mu_buf;
	std::vector<double> xi_s_0_buf;
	std::vector<double> xi_s_2_buf;
	std::vector<double> xi_s_4_buf;
	std::vector<std::vector<double> > xi_wedge_buf;
private:						// Function
	void set_s_mu_buf( double s_min, double s_max,
					   int s_bin_num, int mu_bin_num );
	double inner_integrand( const double & r_sigma,
							const double & r_pi,
							const double & y );
	double inner_integration( const double & s,
							  const double & mu_s );
	double outer_integration( int order, const double & s );
	double wedge_integration( const double & mu_min,
							  const double & mu_max,
							  const double & s );	
public:
	void gen_xi_s(  );
	
	////////// Output //////////
private:
	std::string out_file_name;
public:
	void output(  );

	////////// Misc func/const //////////
private:						// Data
	static const double pi = 3.14159265358979323846;
	static const double nearly0 = 1e-8;
private:						// Function
	integral intg;
	prog_bar pg;
	double legendre( const int & order, const double & x );
};

#endif

