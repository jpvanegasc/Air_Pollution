DEBUG = -g -Wall -fsanitize=address -fsanitize=leak -fsanitize=undefined
VAL = --track-origins=yes --leak-check=full
CACH = --tool=cachegrind
F_STD = -std=f95
OLD = $< -o $@
all : a.out

a.out : main.f90 LB_D3Q19.mod Constants.mod
	gfortran $(F_STD) main.f90 LB_D3Q19.f90 Constants.f90
	time ./$@

LB_D3Q19.mod : LB_D3Q19.f90
	gfortran -c $<

Constants.mod : Constants.f90
	gfortran -c $<

.PHONY : debug
debug : main.f90 LB_D3Q19.mod Constants.mod
	gfortran $(F_STD) $(DEBUG) main.f90 LB_D3Q19.f90 Constants.f90
	./a.out

.PHONY : valgrind
valgrind : main.f90 LB_D3Q19.mod Constants.mod
	g++ $(F_STD) -g main.f90 LB_D3Q19.f90 Constants.f90
	valgrind $(VAL) ./a.out

.PHONY : cachegrind
cachegrind : main.f90 LB_D3Q19.mod Constants.mod
	g++ $(F_STD) -g main.f90 LB_D3Q19.f90 Constants.f90
	valgrind $(CACH) ./a.out

.PHONY : gprof
gprof : main.f90 LB_D3Q19.mod Constants.mod
	g++ $(F_STD) -Wall -pg main.f90 LB_D3Q19.f90 Constants.f90
	./a.out
	gprof a.out gmon.out > analysis.txt
	most analysis.txt

.PHONY : perf
perf : main.f90 LB_D3Q19.mod Constants.mod
	g++ $(F_STD) -Wall -pg main.f90 LB_D3Q19.f90 Constants.f90
	perf record ./a.out ; perf report

.PHONY : clean
clean :
	rm -f *.o *.txt *.out* debug main *.data

.PHONY : plot
plot : data.txt
	gnuplot plot.gp
