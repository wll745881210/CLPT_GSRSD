#ifndef PAIR_V_H
#define PAIR_V_H

#include "corr_func.h"

#include <vector>
#include <string>

class pair_v : public corr_func
{
private:			// Overrider
    void kernel( const double & r, const vec3 & y );
    void post_proc(  );
private:			// Function
    double v12_L( const double & r );
};

#endif

