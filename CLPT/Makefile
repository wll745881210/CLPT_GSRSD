CC = g++-4.8
Option = -O3 -Wall -fopenmp -std=c++0x 
# c++0x or c++11 is for std::to_string on gcc only
# If you are using Clang/llvm, you can remove it

Objs = main.o \
	integral.o lu_decomp.o \
	k_depend_funcs.o \
	q_depend_funcs_single.o q_depend_funcs.o \
	corr_func.o pair_xi.o pair_v.o pair_s.o \
	input.o save_load.o

Out_file = clpt.exe

main : $(Objs)
	$(CC) $(Option) $(Objs) -o $(Out_file)

main.o : main.cpp
	$(CC) $(Option) -c -o main.o main.cpp

corr_func.o : corr_func.cpp corr_func.h
	$(CC) $(Option) -c -o corr_func.o corr_func.cpp

pair_xi.o : pair_xi.cpp pair_xi.h
	$(CC) $(Option) -c -o pair_xi.o pair_xi.cpp

pair_v.o : pair_v.cpp pair_v.h
	$(CC) $(Option) -c -o pair_v.o pair_v.cpp

pair_s.o : pair_s.cpp pair_s.h
	$(CC) $(Option) -c -o pair_s.o pair_s.cpp

integral.o : integral.cpp integral.h
	$(CC) $(Option) -c -o integral.o integral.cpp

lu_decomp.o : lu_decomp.cpp lu_decomp.h
	$(CC) $(Option) -c -o lu_decomp.o lu_decomp.cpp

q_depend_funcs_single.o : \
	q_depend_funcs_single.cpp q_depend_funcs_single.h
	$(CC) $(Option) -c -o q_depend_funcs_single.o \
	q_depend_funcs_single.cpp 

q_depend_funcs.o : q_depend_funcs.cpp q_depend_funcs.h
	$(CC) $(Option) -c -o q_depend_funcs.o \
	q_depend_funcs.cpp 

k_depend_funcs.o : k_depend_funcs.cpp k_depend_funcs.h
	$(CC) $(Option) -c -o k_depend_funcs.o \
	k_depend_funcs.cpp

input.o : input.cpp input.h
	$(CC) $(Option) -c -o input.o input.cpp

save_load.o : save_load.cpp save_load.h
	$(CC) $(Option) -c -o save_load.o save_load.cpp

clean:
	rm *.o *~ clpt.exe
