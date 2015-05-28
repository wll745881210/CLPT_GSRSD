#ifndef PAIR_XI_H
#define PAIR_XI_H

#include "corr_func.h"

class pair_xi : public corr_func
{
private:
    void kernel( const double r, const vec3 & y );
};

#endif
