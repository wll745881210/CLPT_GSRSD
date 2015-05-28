#ifndef LU_DECOMP_H
#define LU_DECOMP_H

class lu_decomp
{
    ////////// Con-/destructor & initializer //////////
public:
    lu_decomp(  );
    ~lu_decomp(  );

    void set_size( int size );

    ////////// Index function //////////
private:
    int l;
public:
    int idx( const int & i, const int & j );

    ////////// LU process //////////
public:
    double lu_inverse( double A[], double A_inv[] );
};

#endif
