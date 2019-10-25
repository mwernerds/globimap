MYBOOSTINC=/home/martin/lib/boost/boost_1_71_0
MYBOOSTLIB=/usr/local/lib
CFLAGS=-std=c++11 -Ofast -fopenmp -march=native -fPIC


globimap.so: globimap.o
	g++ -shared -o globimap.so globimap.o `pkg-config --libs python` -lboost_python27 -lboost_numpy27 -lboost_system -lgomp
globimap.o: globimap.hpp globimap.cpp
	g++ $(CFLAGS) -o globimap.o -c globimap.cpp   `pkg-config --cflags python`

clean:
	rm -f globimap.o globimap.so 

