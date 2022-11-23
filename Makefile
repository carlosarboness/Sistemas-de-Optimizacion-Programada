all: comp exc
comp_check:
	g++ -O2 -Wall -std=c++17 check.cc

comp_exhv:
	g++ -O2 -Wall -std=c++17 exh_vectors.cc

comp_exh: 
	g++ -O2 -Wall -std=c++17 exh.cc

comp_exh2: 
	g++ -O2 -Wall -std=c++17 exh2.cc

comp_exhlb: 
	g++ -O2 -Wall -std=c++17 exh_lb.cc

comp_greedy:
	g++ -O2 -Wall -std=c++17 greedy.cc

comp_cal:
	g++ -O2 -Wall -std=c++17 calculate_penalizations.cc

check:
	./a.out ./med/med-1.txt ./sol_med_exhlb/sol_med-1.txt
	./a.out ./med/med-2.txt ./sol_med_exhlb/sol_med-2.txt
	./a.out ./med/med-3.txt ./sol_med_exhlb/sol_med-3.txt
	./a.out ./med/med-4.txt ./sol_med_exhlb/sol_med-4.txt
	./a.out ./med/med-5.txt ./sol_med_exhlb/sol_med-5.txt
	./a.out ./med/med-6.txt ./sol_med_exhlb/sol_med-6.txt
	./a.out ./med/med-7.txt ./sol_med_exhlb/sol_med-7.txt
	./a.out ./med/med-8.txt ./sol_med_exhlb/sol_med-8.txt
	./a.out ./med/med-9.txt ./sol_med_exhlb/sol_med-9.txt
	./a.out ./med/med-10.txt ./sol_med_exhlb/sol_med-10.txt

exh: 
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
	./a.out ./med/med-1.txt ./sol_med/sol_med-1.txt
	./a.out ./med/med-2.txt ./sol_med/sol_med-2.txt
	./a.out ./med/med-3.txt ./sol_med/sol_med-3.txt
	./a.out ./med/med-4.txt ./sol_med/sol_med-4.txt
	./a.out ./med/med-5.txt ./sol_med/sol_med-5.txt
	./a.out ./med/med-6.txt ./sol_med/sol_med-6.txt
	./a.out ./med/med-7.txt ./sol_med/sol_med-7.txt
	./a.out ./med/med-8.txt ./sol_med/sol_med-8.txt
	./a.out ./med/med-9.txt ./sol_med/sol_med-9.txt
	./a.out ./med/med-10.txt ./sol_med/sol_med-10.txt

cal:
	./a.out ./cal_pen/easy-1.txt

greedy:
	./a.out ./easy/easy-1.txt ./sol_easy_greedy/sol_easy-1.txt
	./a.out ./easy/easy-2.txt ./sol_easy_greedy/sol_easy-2.txt
	./a.out ./easy/easy-3.txt ./sol_easy_greedy/sol_easy-3.txt
	./a.out ./easy/easy-4.txt ./sol_easy_greedy/sol_easy-4.txt
	./a.out ./easy/easy-5.txt ./sol_easy_greedy/sol_easy-5.txt
	./a.out ./easy/easy-6.txt ./sol_easy_greedy/sol_easy-6.txt
	./a.out ./easy/easy-7.txt ./sol_easy_greedy/sol_easy-7.txt
	./a.out ./easy/easy-8.txt ./sol_easy_greedy/sol_easy-8.txt
	./a.out ./easy/easy-9.txt ./sol_easy_greedy/sol_easy-9.txt
	./a.out ./easy/easy-10.txt ./sol_easy_greedy/sol_easy-10.txt	
