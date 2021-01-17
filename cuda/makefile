# Useful flags
VAL = --track-origins=yes --leak-check=full
CACH = --tool=cachegrind
C_STD = -std=c++11

# Dependecies
DEPS = LB_D2Q9.cu

all : main.o

main.o : main.cu $(DEPS) *.h
	nvcc $(C_STD) -O2 $< $(DEPS) -o $@
	time ./$@

.PHONY : valgrind
valgrind : main.cu
	nvcc $(C_STD) -g $< $(DEPS) 
	valgrind $(VAL) ./a.out

.PHONY : cachegrind
cachegrind : main.cu
	nvcc $(C_STD) -g $< $(DEPS) 
	valgrind $(CACH) ./a.out

.PHONY : gprof
gprof : main.cu
	nvcc $(C_STD) -Wall -pg $< $(DEPS) 
	./a.out
	gprof a.out gmon.out > analysis.txt
	most analysis.txt

.PHONY : perf
perf : main.cu
	nvcc $(C_STD) -Wall -pg $< $(DEPS) 
	perf record ./a.out ; perf report

.PHONY : clean
clean :
	rm -f *.o *.txt *.out* debug main *.data

