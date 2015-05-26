CC = g++-4.8
Option = -O3 -Wall -fopenmp

# Objs = main.o corr_func.o pair_v.o pair_s.o integral.o \
# 	lu_decomp.o q_depend_funcs.o k_depend_funcs.o \
# 	input.o prog_bar.o

Objs = main.o integral.o q_depend_funcs_single.o \
	lu_decomp.o q_depend_funcs.o k_depend_funcs.o \
	input.o prog_bar.o

Out_file = clpt.exe

main : $(Objs)
	$(CC) $(Option) $(Objs) -o $(Out_file)

main.o : main.cpp
	$(CC) $(Option) -c -o main.o main.cpp

# corr_func.o : corr_func.cpp corr_func.h
# 	$(CC) $(Option) -c -o corr_func.o corr_func.cpp

# pair_v.o : pair_v.cpp pair_v.h
# 	$(CC) $(Option) -c -o pair_v.o pair_v.cpp

# pair_s.o : pair_s.cpp pair_s.h
# 	$(CC) $(Option) -c -o pair_s.o pair_s.cpp

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

prog_bar.o : prog_bar.cpp prog_bar.h
	$(CC) $(Option) -c -o prog_bar.o prog_bar.cpp

clean:
	rm *.o *~ clpt.exe
