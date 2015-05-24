#include "lu_decomp.h"

#include <cmath>

////////////////////////////////////////////////////////////
// Constructor, desctructor and initializer

lu_decomp::lu_decomp(  )
{
    l = 3;
}

lu_decomp::~lu_decomp(  )
{
	
}

void lu_decomp::set_size( int size )
{
    this->l = size;
    return;
}

////////////////////////////////////////////////////////////
// Index function

int lu_decomp::idx( const int & i, const int & j )
{
    return i * l + j;
}

double lu_decomp::lu_inverse( double A     [  ],
                              double A_inv [  ] )
{
    double det( 1. );

    for( int j = 0; j < l; ++ j )
    {
	for( int i = 0; i < j + 1; ++ i )
	    for( int k = 0; k < i; ++ k )
		A[idx( i, j )] -=
		    A[ idx( i, k ) ] * A[ idx( k, j ) ];
	for( int i = j + 1; i < l; ++ i )
	{
	    for( int k = 0; k < j; ++ k )
		A[idx( i, j )] -=
		    A[ idx( i, k ) ] * A[ idx( k, j ) ];
	    A[ idx( i, j ) ] /= A[ idx( j, j ) ];
	}
	det *= A[idx( j, j )];
    }

    for( int i = 0; i < l; ++ i )
	for( int j = 0; j < l; ++ j )
	    A_inv[idx( i, j )] = ( i==j ? 1. : 0. );

    for( int k = 0; k < l; ++ k )
    {
	for( int i = 1; i < l; ++ i )
	    for( int j = 0; j < i; ++ j )
		A_inv[idx( i, k )] -= A[idx( i, j )]
		    * A_inv[idx( j, k )];
	
	for( int i = l - 1; i > -1; -- i )
	{		
	    for( int j = i + 1; j < l; ++ j )
		A_inv[idx( i, k )] -= A[idx( i, j )]
		    * A_inv[idx( j, k )];
	    A_inv[idx( i, k )] /= A[idx( i, i )];
	}
    }

    return det;
}



