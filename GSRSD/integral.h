#ifndef INTEGRAL_H
#define INTEGRAL_H

class integral
{
	////////// Con-/destructor & initializer//////////
public:
	integral(  );
	~integral(  );
	void clear(  );

	////////// Data buffer //////////
private:
	int counter;
	double x_buf[ 2 ];
	double y_buf[ 2 ];
	double intg_res;

	////////// Read data //////////
public:
	void read( const double & x, const double & y );

	////////// Feed back //////////
public:
	double result(  );

	////////// Gauss-Legendre //////////
public:
	void gl_clear(  );
	double gl_xi( const int & i );
	void gl_read( const double & kernel );
	double gl_result(  );
	static const int gl_num;
private:
	int i_current;
	double gl_intg_res;
	static const double x[  ];
    static const double w[  ];


};

#endif

