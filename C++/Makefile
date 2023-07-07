OPTIONS = -std=c++17 -O3 -DNDEBUG

binaries = tests demo demo_ms convergence_test boundary_test

all: ${binaries}

tests: tests.cpp
	g++ ${OPTIONS} -o $@ $^

demo: demo.cpp
	g++ ${OPTIONS} -o $@ $^

demo_ms: demo_ms.cpp
	g++ ${OPTIONS} -o $@ $^

convergence_test: convergence_test.cpp
	g++ ${OPTIONS} -o $@ $^

boundary_test: boundary_test.cpp
	g++ ${OPTIONS} -o $@ $^	


clean:
	rm -f ${binaries}
