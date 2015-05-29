CC = g++-4.8
Option = -O3 -Wall -fopenmp -std=c++11

Objs = main.o integral.o xi_stream.o interp1.o input.o

Out_file = xi_stream.exe

main : $(Objs)
	$(CC) $(Option) $(Objs) -o $(Out_file)

main.o : main.cpp
	$(CC) $(Option) -c -o main.o main.cpp

integral.o : integral.cpp integral.h
	$(CC) $(Option) -c -o integral.o integral.cpp

xi_stream.o : xi_stream.cpp xi_stream.h
	$(CC) $(Option) -c -o xi_stream.o xi_stream.cpp

interp1.o : interp1.cpp interp1.h
	$(CC) $(Option) -c -o interp1.o interp1.cpp

input.o : input.cpp input.h
	$(CC) $(Option) -c -o input.o input.cpp

clean:
	rm *.o *~ *.exe

