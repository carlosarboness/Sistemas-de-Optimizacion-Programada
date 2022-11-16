all: comp exc
comp: 
	g++ -O2 -std=c++17 exh.cc
exc: 
	./a.out in.txt out.txt