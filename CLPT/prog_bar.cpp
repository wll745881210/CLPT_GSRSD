#include <iostream>
#include "prog_bar.h"

////////////////////////////////////////////////////////////
// Constructor and destructor

prog_bar::prog_bar(  )
{
	
}

prog_bar::~prog_bar(  )
{
	
}

void prog_bar::init( const unsigned & total )
{
	this->total = total - 1;
	this->lock  = false;
	this->last_k = 0;
	return;
}

void prog_bar::show( const unsigned & current )
{
	const unsigned k = 10 * current / total;

	if( lock && k > 0 )
		lock = false;
	else if( k != last_k && ( !lock ) )
	{
		std::cout << "=>" << k * 10 << '%';
		std::cout.flush(  );
		lock = true;
		if( k == 10 )
			std::cout << std::endl;
	}
	last_k = k;
	
	return;
}
