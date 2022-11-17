all: comp exc
comp: 
	g++ -O2 -std=c++17 exh.cc
exc: 
	./a.out ./easy/easy-1.txt ./sol_easy/sol_easy-1.txt
	./a.out ./easy/easy-2.txt ./sol_easy/sol_easy-2.txt
	./a.out ./easy/easy-3.txt ./sol_easy/sol_easy-3.txt
	./a.out ./easy/easy-4.txt ./sol_easy/sol_easy-4.txt
	./a.out ./easy/easy-5.txt ./sol_easy/sol_easy-5.txt
	./a.out ./easy/easy-6.txt ./sol_easy/sol_easy-6.txt
	./a.out ./easy/easy-7.txt ./sol_easy/sol_easy-7.txt
	./a.out ./easy/easy-8.txt ./sol_easy/sol_easy-8.txt
	./a.out ./easy/easy-9.txt ./sol_easy/sol_easy-9.txt
	./a.out ./easy/easy-10.txt ./sol_easy/sol_easy-10.txt

cal:
	g++ -O2 -std=c++17 calculate_penalizations.cc
	./a.out ./cal_pen/easy-1.txt
