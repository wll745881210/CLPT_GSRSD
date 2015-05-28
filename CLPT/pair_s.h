#ifndef PAIR_S_H
#define PAIR_S_H

#include "corr_func.h"

class pair_s : public corr_func
{
public:				// Override the constructor
    pair_s(  );
private:			// Override
    void kernel( const double & r, const vec3 & y,
                 double * bias_comp );
    void post_proc(  );
};

#endif
