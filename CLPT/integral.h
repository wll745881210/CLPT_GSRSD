#ifndef INTEGRAL_H
#define INTEGRAL_H

#include <vector>

class integral
{
    ////////// Con-/destructor & initializer//////////
public:
    integral(  );
    ~integral(  );
    void clear(  );

    ////////// Data buffer //////////
private:
    std::vector<double> x_buf, y_buf;

    ////////// Read data //////////
public:
    void read( const double & x, const double & y );

    ////////// Feed back //////////
public:
    double result(  );

    ////////// Gauss-Legendre //////////
public:
    void   gl_clear(  );
    double gl_xi( const int & i );
    void   gl_read( const int & i, const double & kernel );
    double gl_result(  );
    static const int gl_num;
private:
    double gl_intg_res;
    static const double x[  ];
    static const double w[  ];
};

#endif

