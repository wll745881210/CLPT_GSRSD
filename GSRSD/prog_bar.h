#ifndef PROG_BAR_H
#define PROG_BAR_H


class prog_bar
{
	////////// Con-/destructor //////////
public:
	prog_bar(  );
	~prog_bar(  );
	
	////////// Status bar //////////
public:
	void init( const unsigned & total );
	void show( const unsigned & current );
private:
	unsigned total;
	bool lock;
	unsigned last_k;
};

#endif
